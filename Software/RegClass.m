classdef RegClass < handle
% REGCLASS Class containing methods to find optimal registration.
%
% Description:
%       This class is used for finding the optimal registration between
%   two images of different modalities. The registration method is based
%   on maximization of mutual information. The SAM image is stored in the
%   ImageClass property 'img_a', the histology image is stored in the
%   ImageClass property 'img_b'. The class searches for the SAM image
%   within the histology image, and therefore it is assumed that the
%   histology image is larger than the SAM image.
%

    properties
        img_a                       % SAM image ImageClass
        img_b                       % Histology image ImageClass
        hist                        % 3 Channel Histology
        joint_pmf                   % joint pmf of img_a and img_b
        joint_entropy               % joint entropy of img_a and img_b
        mutual_info                 % mutual information btwn img_a, img_b
        bin_sz                      % bin size used in entropy calculation
        normalized_mutual_info      % NMI between img_a, img_b
    end

    methods
        function obj = RegClass(a, b, a_d, b_d)
            % REGCLASS Constructor for the registration class.
            % 'a' is SAM image, 'b' is histology, bin_sz is
            % used for entropy calculations. After these properties are
            % initalized, we calculate the joint pmf, joint entropy, mutual
            % information, and normalized mutual information. All
            % properties are saved within.
            %
            
            if nargin == 4
                %% preprocess image data
                % make image size an integer multiple of 16 (because we 
                % downsample by 16)
                % NOTE: We only register grayscale images!
                trim = floor(size(b(:, :, 1))/16);
                obj.hist = b(1:trim(1)*16, 1:trim(2)*16, :);
                b = b(1:trim(1)*16, 1:trim(2)*16, 1); % only grayscale
                trim = floor(size(a)/16);
                a = a(1:trim(1)*16, 1:trim(2)*16);
                
                %% initalize properties 
                obj.bin_sz = 1;
                obj.img_a = ImageClass(a, obj.bin_sz, a_d);
                obj.img_b = ImageClass(b, obj.bin_sz, b_d, ...
                                       size(obj.img_a.img_full));
                obj.FormJointPMF();
                obj.FormJointEntropy();
                obj.FormMutualInfo();
                obj.FormNormalizedMutualInfo();
            end
        end
        function obj = UpdateReg(obj)
            % UPDATEREG Update all properties of the object.
            % Updates image properties, then updates all of the
            % registration properties.
            %
            
            obj.img_a.UpdateImg();
            obj.img_b.UpdateImg();
            obj.FormJointPMF();
            obj.FormJointEntropy();
            obj.FormMutualInfo();
            obj.FormNormalizedMutualInfo();
        end
        function val = HandleOpt(obj, x_T, y_T, t_T)
            % HANDLEOPT This function is a wrapper for fminsearch.
            % For quicker operation, two modes of transformation are
            % available: a complete transform with a translation and
            % rotation, and a simplified transform with only a translation.
            % This reduces the complexity of our optimization scheme
            % (searching over a two parameter range is faster than
            % searching over a three parameter range). This simplification
            % exploits the prior that the angle between the histology and
            % the SAM images is quite small.
            %
            
            % if an angle is passed in, do the whole transformation (theta)
            if nargin == 4
                obj.img_b.TransformImage([x_T, y_T, t_T]);
            % otherwise just translate
            else
                obj.img_b.TranslateImage([round(x_T), round(y_T)]);
            end
            
            % update
            obj.UpdateReg();
            
            % return similarity metric (nmi) to maximize
            val = obj.normalized_mutual_info;
        end
        function obj = FormJointPMF(obj)
            % FORMJOINTPMF Calculates the joint probability mass function
            % of the two images.
            % 
            
            obj.joint_pmf = histcounts2(obj.img_a.img, ...
                obj.img_b.img, ...
                [256 / obj.bin_sz, 256 / obj.bin_sz], ...
                'Normalization', 'probability');
        end
        function trans_pmf = ShowJointPMF(obj)
            % SHOWJOINTPMF This function is used to display the joint PMF.
            % It transforms the joint PMF simply to make the visualization
            % of it more clear.
            %
            
            const = 255 / max(obj.joint_pmf(:));
            trans_pmf = const * obj.joint_pmf;
            imshow(trans_pmf);
        end
        function obj = FormJointEntropy(obj)
            % FORMJOINTENTROPY Simply calculates the joint entropy between
            % two images.
            % This implementation is inefficient and can be vectorized.
            %
            
            jnt_entr = 0;
            for i=1:256 / obj.bin_sz
                for j=1:256 / obj.bin_sz
                    temp = log2(obj.joint_pmf(i, j));
                    temp(isinf(temp)) = 0;
                    temp = - obj.joint_pmf(i, j) * temp;
                    jnt_entr = jnt_entr + temp;
                end
            end
            obj.joint_entropy = jnt_entr;
        end
        function obj = FormMutualInfo(obj)
            % FORMMUTUALINFO Calculates the mutual information between the 
            % two images.
            %
            
            obj.mutual_info = obj.img_a.entropy + obj.img_b.entropy ...
                - obj.joint_entropy;
        end
        function obj = FormNormalizedMutualInfo(obj)
            % FORMNORMALIZEDMUTUALINFO Calculates the NMI between the two 
            % images.
            %
            
            obj.normalized_mutual_info = (obj.img_a.entropy + ...
                obj.img_b.entropy) / obj.joint_entropy;
        end
        
        function start_point = FindStart(obj, num_s)
            %START_POINT = FIND_START(NUM_S) finds an initial starting
            % point on a downscaled histology image.
            % Description:
            % A high level summary is as follows:
            %  - Downsample images.
            %  - Change pixel size of histology image to match SAM image.
            %  - Uniformly generate random samples throughout histology image.
            %  - Filter the random sample points.
            %  - Run a simplex optimization at each of the remaining samples.
            %  - Choose the resulting optima with the largest NMI for start.
            %  - Convert the optimal transform to original size and rate.
            
            % temporarily push downsampled image at original rate
            ds_img = obj.img_b.og_img;
            og_dx = obj.img_b.dx;
            og_dy = obj.img_b.dy;
            
            % ratechange downsampled image to SAM rate, & find start point
            obj.img_b.RateChange([obj.img_a.dx, obj.img_a.dy]);
            
            %% random sampling
            % boundaries on samples
            x_max = size(obj.img_b.og_img, 2);
            y_max = size(obj.img_b.og_img, 1);

            % draw samples
            x_samps = randi(x_max, [num_s, 1]);
            y_samps = randi(y_max, [num_s, 1]);

            % filter sampled points based on entropy
            entr_thresh = 5;
            y_center = ceil(size(obj.img_b.og_img, 1)/2);
            x_center = ceil(size(obj.img_b.og_img, 2)/2);
            x_s = x_center;
            y_s = y_center;
            for i=1:num_s
                obj.HandleOpt(x_center - x_samps(i), ...
                              y_center - y_samps(i));
                if (obj.img_b.entropy > entr_thresh)
                    % if it's a good sample, save it
                    x_s = [x_samps(i), x_s];
                    y_s = [y_samps(i), y_s];
                end
            end


            %% evaluate NMI at sampled points
            fun = @(x)(2 - obj.HandleOpt(x(1), x(2), x(3)));
            opt_loc = cell(size(x_s, 2), 1);
            nmi = zeros(size(x_s, 2), 1);
            num = size(x_s, 2);
            options = optimset('MaxIter', 15, 'TolFun', 1e-3, 'Display', 'off');
            for i=1:size(x_s, 2)
                opt_T = fminsearch(fun, [x_center - x_s(i), ...
                                         y_center - y_s(i), ...
                                         -2], ...
                                         options);
                opt_loc{i} = opt_T;
                nmi(i) = obj.HandleOpt(opt_T(1), opt_T(2), opt_T(3));
                if mod(i, 5) == 0
                    fprintf('finding starting point... %d of %d\n', i, num);
                end
            end

            % find norm between the two best estimates
            [sorted_nmi, ind] = sort(nmi, 'descend');
            dist = norm([opt_loc{ind(1)}(1), ...
                         opt_loc{ind(1)}(2)] - ...
                         [opt_loc{ind(2)}(1), ...
                         opt_loc{ind(2)}(2)]);

            % warn user that they may have to rerun the script
            if (dist > 30) 
                fprintf('large norm between two best estimates... resample! \n')
            end

            % the best starting point
            start_point = [opt_loc{ind(1)}(1), ...
                           opt_loc{ind(1)}(2), ...
                           opt_loc{ind(1)}(3)];

            % update image to best location
            obj.HandleOpt(start_point(1), start_point(2), start_point(3));

            % starting point for full size image (at matched sampling rate)
            disp([start_point(1), start_point(2)])

            % pop downsampled image to original rate
            obj.img_b.og_img = ds_img;
            obj.img_b.dx = og_dx;
            obj.img_b.dy = og_dy;
            
            % convert starting points back to old sampling rate
            start_point(1) = start_point(1) * obj.img_a.dx/obj.img_b.dx;
            start_point(2) = start_point(2) * obj.img_a.dy/obj.img_b.dy;

            % back to full size (and original rate)
            start_point(1:2) = ceil(16*start_point(1:2));

        end
        
        function hist = Register(obj)
            %HIST = REGISTER registers hist image with SAM amplitude image.
            %   
            % This function attempts to register images of different modalities.
            % It is not guaranteed to find the optimal registration, so it necessary
            % to manually inspect the coregistration after completion. It may be
            % necessary to rerun this function if the registration is not
            % satisfactoy.

            %% find starting point
            % number of starting points
            bb = size(obj.img_b.og_img, 1) * size(obj.img_b.og_img, 2);
            s_num = ceil(log(0.05)/log((bb - 65*65)/bb));
            fprintf('evaluating %d starting samples... \n', s_num);
            start_point = obj.FindStart(s_num);
            
            %% trim histology according to start point (speed up interpolation)
            % we add two borders... 
            % the first border makes up for the optimal angle being set to 0, and 
            % the second border makes up for a possible rotation of 45 degrees

            % for angle
            H = sqrt(start_point(1)^2 + start_point(2)^2);
            y_c = size(obj.img_b.og_img_full, 1)/2;
            y_lo = y_c - (start_point(2) + ...
                          abs(round(H*sind(start_point(3))))) - ...
                   floor(size(obj.img_a.img_full, 1)/2);
            y_hi = y_c - (start_point(2) - ...
                          abs(round(H*sind(start_point(3))))) + ...
                   floor(size(obj.img_a.img_full, 1)/2);
            x_c = ceil(size(obj.img_b.og_img_full, 2)/2);
            x_lo = x_c - (start_point(1) + ...
                          abs(round(H*(1 - cosd(start_point(3)))))) - ...
                   floor(size(obj.img_a.img_full, 2)/2);
            x_hi = x_c - (start_point(1) - ...
                          abs(round(H*(1 - cosd(start_point(3)))))) + ...
                   floor(size(obj.img_a.img_full, 2)/2);

            % for rotation
            y_sz = ceil(0.3*size(y_lo:y_hi, 2));
            x_sz = ceil(0.3*size(x_lo:x_hi, 2));

            % update downsampled histology with trimmed image
            obj.img_b.og_img = obj.img_b.og_img_full( ...
                                          max(y_lo - y_sz, 1) : ...
                                          min(y_hi + y_sz, 2*y_c), ...
                                          max(x_lo - x_sz, 1) : ...
                                          min(x_hi + x_sz, 2*x_c), ...
                                          :);
            trim = floor(size(obj.img_b.og_img(:, :, 1))/16);
            obj.img_b.og_img = obj.img_b.og_img(1:trim(1)*16, ...
                                                1:trim(2)*16, :);

            % 3 channel hist
            obj.hist = obj.hist( ...
                                max(y_lo - y_sz, 1) : ...
                                min(y_hi + y_sz, 2*y_c), ...
                                max(x_lo - x_sz, 1) : ...
                                min(x_hi + x_sz, 2*x_c), ...
                                :);
            trim = floor(size(obj.hist(:, :, 1))/16);
            obj.hist = obj.hist(1:trim(1)*16, ...
                                1:trim(2)*16, :);
            
            % rate change entire histology
            hist = obj.RateChange3Channel(obj.hist, ...
                                          [obj.img_b.dx, obj.img_b.dy], ...
                                          [obj.img_a.dx, obj.img_a.dy]);

            % update downsampled SAM with original image
            obj.img_a.img = obj.img_a.img_full;
            
            % change sampling rate of histology
            obj.img_b.RateChange([obj.img_a.dx, obj.img_a.dy]);
            
            % reset size of obj.img_b.img
            obj.img_b.img = zeros(size(obj.img_a.img));
            
            
            % update properties
            obj.UpdateReg();
            
            %% final optimization
            % function handle for fminsearch
            fun = @(x)(2 - obj.HandleOpt(x(1), x(2), x(3)));

            % need to initialize with nonzero values (fminsearch calculates
            % step size based on a percentage of the initial values)
            opt_T = fminsearch(fun, [100, 100, sign(start_point(3))*4])

            % final registration
            obj.HandleOpt(opt_T(1), opt_T(2), opt_T(3));

            % transform entire histology... should make a class for this
            hist = obj.Transform3Channel(hist, opt_T, size(obj.img_a.img));
        end
        
        function n_img = RateChange3Channel(o_img, o_dx_dy, n_dx_dy)

            % get size of image at original rate
            o_y_max = size(o_img, 1);
            o_x_max = size(o_img, 2);

            % find size (in microns) of image
            y_max = o_y_max * o_dx_dy(2);
            x_max = o_x_max * o_dx_dy(1);

            % original image is defined on this grid
            [X, Y] = meshgrid([0:o_dx_dy(1):x_max - o_dx_dy(1)], ...
                              [0:o_dx_dy(2):y_max - o_dx_dy(2)]);

            % new image will be defined on this grid
            [Xq, Yq] = meshgrid([0:n_dx_dy(1):round(x_max)], ...
                                [0:n_dx_dy(2):round(y_max)]);

            % rate changed image
            n_img = zeros(size(Xq, 1), size(Xq, 2), 3);
            for i=1:3
                n_img(:, :, i) = interp2(X, Y, single(o_img(:, :, i)), Xq, Yq);
            end
            n_img = uint8(n_img); 
        end
        
        function img = Transform3Channel(img, params, n_img_sz)

            % get transformation parameters into correct data structure
            params = num2cell(params);
            [T_x, T_y, theta] = params{:};

            % transformation matrix
            A = [cosd(theta), -sind(theta), T_x; ...
                 sind(theta), cosd(theta), T_y; ...
                 0, 0, 1];

            % form coordinate matrix (center at <0, 0>)
            x_sz = size(img, 2);
            y_sz = size(img, 1);
            [x, y] = meshgrid([-floor(x_sz/2):ceil(x_sz/2 - 1)], ... 
                              [-floor(y_sz/2):ceil(y_sz/2 - 1)]);

            % coordinates for the transformed image to be defined at
            num = numel(x);
            new_loc =  [reshape(x, [1, num]); ...
                        reshape(y, [1, num]); ...
                        ones(1, num)];

            % inverse mapping function to find original pixel locations
            loc = A\new_loc; % inv(A)*b
            u = reshape(loc(1, :), size(img(:, :, 1)));
            v = reshape(loc(2, :), size(img(:, :, 1)));

            % warp image to fit onto specified grid
            T_img = zeros(size(u, 1), size(u, 2), 3);
            for i=1:3
                T_img(:, :, i) = interp2(x, y, single(img(:,:, i)), u, v, 'linear');
            end

            % trim back to smaller size (at center!)
            center_y = round(size(img, 1)/2);
            center_x = round(size(img, 2)/2);
            row_lo = center_y - floor(n_img_sz(1)/2);
            row_hi = row_lo + n_img_sz(1) - 1;
            col_lo = center_x - floor(n_img_sz(2)/2);
            col_hi = col_lo + n_img_sz(2) - 1;
            img = uint8(T_img(row_lo:row_hi, col_lo:col_hi, :));
        end
        
    end
end


