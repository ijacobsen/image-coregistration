classdef ImageClass < handle
% IMAGECLASS A class to store and handle images and associated metrics
% that are used in image registration.
%
% Description:
%       This class stores images as arrays, and calculates useful metrics
%   about the images. These metrics are used to perform image registration
%   with the RegClass class.

    properties
        img;            % image array
        img_full;       % full scale image array
        pmf;            % probability mass function of the image
        entropy;        % entropy of the image
        bin_sz;         % bin size, used for PMF
        og_img;         % original (full size) image array
        og_img_full;    % fullscale original image array
        dx;             % pixel size in x direction
        dy;             % pixel size in y direction
    end
    methods
        function obj = ImageClass(data_in, bin_sz, dx_dy, img_sz)
            % IMAGECLASS Constructor for the ImageClass.
            % 'data_in' is either a filename, or an array containing the
            % image data. 'bin_sz' is the bin size used for the PMF
            % calculation. 'dx_dy' is a vector containing [dx, dy]. 
            % 'img_sz' is an optional input. It is used to specify the size
            % of the SAM image when the currect object is used for a
            % histology image.
            %
            
            if nargin > 1
                if ischar(data_in)
                    obj.og_img_full = imread(data_in);
                elseif isnumeric(data_in)
                    obj.og_img_full = data_in;
                else
                    error('must pass in filename as string or numeric array')
                end 
                if nargin == 4
                    obj.img_full = obj.og_img_full(1:img_sz(1), 1:img_sz(2));
                else
                    obj.img_full = obj.og_img_full;
                end
                obj.img = downsample(downsample(obj.img_full, 16)', 16)';
                obj.og_img = downsample(downsample(obj.og_img_full, 16)', 16)';
                obj.dx = dx_dy(1);
                obj.dy = dx_dy(2);
                obj.bin_sz = bin_sz;
                obj.UpdateImg();
            end
        end
        function obj = RateChange(obj, n_dx_dy)
            % RATECHANGE Changes the sampling rate of an image.
            % 'n_dx_dy' is a vector containing the new sampling rate, of
            % the form [dx, dy]. The function uses linear interpolation in
            % the rate changing.
            %
            
            % get size of image at original rate
            [o_y_max, o_x_max] = size(obj.og_img);

            % find size (in microns) of image
            y_max = o_y_max * obj.dy;
            x_max = o_x_max * obj.dx;

            % original image is defined on this grid
            [X, Y] = meshgrid([0:obj.dx:x_max - obj.dx], ...
                              [0:obj.dy:y_max - obj.dy]);

            % new image will be defined on this grid
            [Xq, Yq] = meshgrid([0:n_dx_dy(1):round(x_max)], ...
                                [0:n_dx_dy(2):round(y_max)]);

            % rate changed image
            n_img = interp2(X, Y, single(obj.og_img), Xq, Yq);
            
            % update image and descriptors
            obj.og_img = uint8(n_img);
            obj.img = obj.og_img(1:size(obj.img, 1), 1:size(obj.img, 2));
            obj.dx = n_dx_dy(1);
            obj.dy = n_dx_dy(2);
            obj.UpdateImg();
            
        end
        function obj = FormPMF(obj)
            % FORMPMF Forms the probability mass function of the image.
            %
            
            obj.pmf = histcounts(obj.img, 256 / obj.bin_sz, ...
                                 'Normalization', 'probability');
        end
        function obj = CalcEntropy(obj)
            % CALCENTOPY Calculates the entropy of the image.
            %
            
            % this can be vectorized
            entr = 0;
            for i=1:256 / obj.bin_sz
                temp=log2(obj.pmf(i));
                temp(isinf(temp)) = 0;
                temp = - obj.pmf(i) * temp;
                entr = entr + temp;
            end
            obj.entropy = entr;
        end
        function obj = UpdateImg(obj)
            % UPDATEIMG Updates the PMF and entropy properties of the
            % image.
            %
            
            % this function should be called after any changes are made
            obj.FormPMF();
            obj.CalcEntropy();
        end
        
        function obj = TranslateImage(obj, params)
            % TRANSLATEIMAGE Translates the image by integer values. 
            % This is only used when the ImageClass object contains the
            % histology image.
            %
            
            % get transformation parameters into correct data structure
            params = num2cell(params);
            [T_x, T_y] = params{:};     

            % translate image in y direction
            temp = circshift(obj.og_img, T_y, 1);

            % translate image in x direction
            temp = circshift(temp, T_x, 2);

            % trim image
            center_y = round(size(obj.og_img, 1)/2);
            center_x = round(size(obj.og_img, 2)/2);
            row_lo = center_y - floor(size(obj.img, 1)/2);
            row_hi = row_lo + size(obj.img, 1) - 1;
            col_lo = center_x - floor(size(obj.img, 2)/2);
            col_hi = col_lo + size(obj.img, 2) - 1;
            obj.img = temp(row_lo:row_hi, col_lo:col_hi);
            
            % update statistics
            obj.UpdateImg();
        end
            

        function obj = TransformImage(obj, params)
            % TRANSFORMIMAGE Rotates and translates the image.
            % This is only used when the ImageClass object contains the
            % histology image.
            %
            
            % get transformation parameters into correct data structure
            params = num2cell(params);
            [T_x, T_y, theta] = params{:};
             
            % transformation matrix
            A = [cosd(theta), -sind(theta), T_x; ...
                 sind(theta), cosd(theta), T_y; ...
                 0, 0, 1];
            
            % form coordinate matrix (center at <0, 0>)
            x_sz = size(obj.og_img, 2);
            y_sz = size(obj.og_img, 1);
            [x, y] = meshgrid([-floor(x_sz/2):ceil(x_sz/2 - 1)], ... 
                              [-floor(y_sz/2):ceil(y_sz/2 - 1)]);
            
            % coordinates for the transformed image to be defined at
            num = numel(x);
            new_loc =  [reshape(x, [1, num]); ...
                        reshape(y, [1, num]); ...
                        ones(1, num)];

            % inverse mapping function to find original pixel locations
            loc = A\new_loc; % inv(A)*b
            u = reshape(loc(1, :), size(obj.og_img));
            v = reshape(loc(2, :), size(obj.og_img));

            % warp image to fit onto specified grid
            T_img = interp2(x, y, double(obj.og_img), u, v, 'linear');
            
            % trim back to smaller size (at center!)
            center_y = round(size(obj.og_img, 1)/2);
            center_x = round(size(obj.og_img, 2)/2);
            row_lo = center_y - floor(size(obj.img, 1)/2);
            row_hi = row_lo + size(obj.img, 1) - 1;
            col_lo = center_x - floor(size(obj.img, 2)/2);
            col_hi = col_lo + size(obj.img, 2) - 1;
            obj.img = uint8(T_img(row_lo:row_hi, col_lo:col_hi));
            
            % update statistics
            obj.UpdateImg();

        end
    end
end
