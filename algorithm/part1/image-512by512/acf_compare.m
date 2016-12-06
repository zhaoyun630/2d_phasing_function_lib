% This piece of code is used to demonstrate autocorrelation function and
% the inverse fourier transform of diffraction pattern.

% Runs very slow, DONT USE THIS CODE

% Written by Yun Zhao
% First created on Aug 25, 2014
% Last editted on Aug 25, 2014

% load the obj
prompt = 'please give the input obj file name: \n';
obj_file_name = input(prompt,'s');
obj = importdata(obj_file_name);

% calculate autocorrelation function(acf) from density
fprintf('calculate autocorrelation function(acf) from density: \n');
[obj_size_a, obj_size_b] = size(obj);
obj_center_a = int16(obj_size_a/2);
obj_center_b = int16(obj_size_b/2);
acf = zeros(2*obj_size_a,2*obj_size_b); % autocorrelation function

%pad zeros around the object in order to simply the calculation for acf
obj_1 = zeros(2*obj_size_a,2*obj_size_b);
for i= 1:obj_size_a
    for j = 1:obj_size_b
        obj_1(i+obj_center_a,j+obj_center_b)= obj(i,j);
    end
end

fprintf(' calculating first quadrant of acf: r0_a >= 0 && r0_b >= 0 \n');
for i = 2*obj_center_a:2*obj_size_a
    for j = 2*obj_center_b:2*obj_size_b
        r0_a = i - 2*obj_center_a;
        r0_b = j - 2*obj_center_b;
        overlap_1 = obj_1(1:(2*obj_size_a-r0_a),1:(2*obj_size_b-r0_b));
        overlap_2 = obj_1(r0_a+1:2*obj_size_a,r0_b+1:2*obj_size_b);
        acf(i,j) = sum(sum(overlap_1.*overlap_2));
    end
end

fprintf('calculating second quadrant of acf: r0_a <= 0 && r0_b >= 0 \n');
for i = 1:2*obj_center_a
    for j = 2*obj_center_b:2*obj_size_b
        r0_a = i - 2*obj_center_a;
        r0_b = j - 2*obj_center_b;
        overlap_1 = obj_1(1:(2*obj_size_a+r0_a),r0_b+1:2*obj_size_b);
        overlap_2 = obj_1(-r0_a+1:2*obj_size_a,1:(2*obj_size_b-r0_b));
        acf(i,j) = sum(sum(overlap_1.*overlap_2));
    end
end

fprintf('calculating third quadrant of acf: r0_a <= 0 && r0_b <= 0 \n');
for i = 1:2*obj_center_a
    for j = 1:2*obj_center_b
        r0_a = i - 2*obj_center_a;
        r0_b = j - 2*obj_center_b;
        overlap_1 = obj_1(1:(2*obj_size_a+r0_a),1:(2*obj_size_b+r0_b));
        overlap_2 = obj_1(-r0_a+1:2*obj_size_a,-r0_b+1:2*obj_size_b);
        acf(i,j) = sum(sum(overlap_1.*overlap_2));
    end
end

fprintf('calculating four quadrant of acf: r0_a >= 0 && r0_b <= 0 \n');
for i = 2*obj_center_a:2*obj_size_a
    for j = 1:2*obj_center_b
        r0_a = i - 2*obj_center_a;
        r0_b = j - 2*obj_center_b;
        overlap_1 = obj_1(1:(2*obj_size_a-r0_a),-r0_b+1:2*obj_size_b);
        overlap_2 = obj_1(r0_a+1:2*obj_size_a,1:(2*obj_size_b+r0_b));
        acf(i,j) = sum(sum(overlap_1.*overlap_2));
    end
end
save('acf_1.mat','acf')



% calculate diffraction pattern.
fprintf('start calculating diffraction pattern: \n');
obj_fft = fft2(obj);
diff_pattern = abs(fftshift(obj_fft)).^2;
diff_pattern_shift = abs(obj_fft).^2;


% calculate acf from diffraction pattern.
fprintf('start calculating acf from diffraction pattern: \n');
acf_shift = abs(fftshift(ifft2(diff_pattern)));
acf_diff = abs(ifft2(diff_pattern_shift));
save('acf_shift.mat','acf_shift');
save('acf_diff.mat','acf_diff');


% for i = 1:2*obj_size_a
%     for j = 1:2*obj_size_b
%         r0_a = i - 2*obj_center_a;
%         r0_b = j - 2*obj_center_b;
%         if r0_a >= 0 && r0_b >= 0
%             overlap_1 = obj_1(1:(2*obj_size_a-r0_a),1:(2*obj_size_b-r0_b));
%             overlap_2 = obj_1(r0_a+1:2*obj_size_a,r0_b+1:2*obj_size_b);
%             acf(i,j) = sum(sum(overlap_1.*overlap_2));
%         elseif r0_a >= 0 && r0_b <= 0
%             overlap_1 = obj_1(1:(2*obj_size_a-r0_a),-r0_b+1:2*obj_size_b);
%             overlap_2 = obj_1(r0_a+1:2*obj_size_a,1:(2*obj_size_b+r0_b));
%             acf(i,j) = sum(sum(overlap_1.*overlap_2));
%         elseif r0_a <= 0 && r0_b <= 0
%             overlap_1 = obj_1(1:(2*obj_size_a+r0_a),1:(2*obj_size_b+r0_b));
%             overlap_2 = obj_1(-r0_a+1:2*obj_size_a,-r0_b+1:2*obj_size_b);
%             acf(i,j) = sum(sum(overlap_1.*overlap_2));
%         elseif r0_a <= 0 && r0_b >= 0
%             overlap_1 = obj_1(1:(2*obj_size_a+r0_a),r0_b+1:2*obj_size_b);
%             overlap_2 = obj_1(-r0_a+1:2*obj_size_a,1:(2*obj_size_b-r0_b));
%             acf(i,j) = sum(sum(overlap_1.*overlap_2));
%         end
%     end
% end        