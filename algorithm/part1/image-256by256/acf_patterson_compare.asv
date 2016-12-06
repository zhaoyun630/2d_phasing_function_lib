% This code is used to calculate Patterson function.
% Then we can compare it with autocorrelation function from definition


prompt = 'Please give an object: \n';
obj_file = input(prompt,'s');
obj = importdata(obj_file);

% Autocorrelation function from definition
acf_def = xcorr2(obj);
save('acf_def.mat','acf_def');

% Now Calculate Patterson function
F_complex = fft2(obj);
F_modules = abs(F_complex);
I = F_modules.^2;
Patterson_function = fftshift(abs(ifft2(I)));
save('Patt_fun.mat','Patterson_function');

figure(1)
imshow(acf_def);

figure(2)
imshow(Patterson_function);


