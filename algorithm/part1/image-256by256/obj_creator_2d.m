% make obj from a image

% first created by Yun Zhao on Aug 25, 2014
% Last editted by Yun Zhao on Aug 25, 2014


prompt = 'provide a image file: \n';
image_file_name = input(prompt,'s');
protein_image = imread(image_file_name);
image_size = size(protein_image);
if image_size(2)==3
    protein_image = rgb2gray(protein_image);
end


prompt = 'provide a support file: \n';
support_file_name = input(prompt,'s');
support = importdata(support_file_name);
support_size = size(support);

fprintf('Now we are going crop the image to the support size \n');
protein_image_crop = zeros(support_size);
center_shift_a = int16((image_size(1)-support_size(1))/2);
center_shift_b = int16((image_size(2)-support_size(1))/2);
for i= 1:support_size(1)
    for j = 1:support_size(2)
        protein_image_crop(i,j) = protein_image(i+center_shift_a,j+center_shift_b)/255;
    end
end

fprintf('Here we apply the support on new cropped image: \n');
fprintf('It is going to be an obj file with zero values more than \n half of its area \n');

obj = support.* protein_image_crop;
save('obj_2d.mat','obj');

fprintf('Here we can check the new image \n')
imshow(obj);
