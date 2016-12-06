% This piece of code is used to calculate fractional volume over 
% sigma level

% First created by Yun Zhao on Sep 8,2014
% Last editted by Yun Zhao on Sep 8,2014

fprintf(' make sure the input obj file is a single unit cell \n');
fprintf(' otherwise, it wouldn''t work properly. \n')
prompts = 'please give the obj file: \n';
obj_file = input(prompts,'s');
obj = importdata(obj_file);
[s_a,s_b,s_c] = size(obj);
N_total = s_a*s_b*s_c;
n_sample = 1000;
sigma_level_max = 10;
volume_frac = zeros(1,n_sample);

for i = 1:n_sample
    sigma_level = (i/n_sample)*sigma_level_max;
    rho_1D_array = reshape(obj,1,s_a*s_b*s_c); % Squeeze 3d array to 1d array
    rho_mean = mean(rho_1D_array);
    rho_sigma = std(rho_1D_array);
    
    support = zeros(s_a,s_b,s_c);
    rho_cutoff = sigma_level*rho_sigma;
    for l=1:s_a
        for m = 1:s_b
            for n=1:s_c
                if obj(l,m,n) < rho_cutoff
                    support(l,m,n) = 0;
                else
                    support(l,m,n) = 1;
                end
            end
        end
    end
    volume_frac(i)= sum(sum(sum(support)))/N_total;
end

x_sigma = linspace(ceil(sigma_level_max)/n_sample,ceil(sigma_level_max),n_sample);


% plt1 = plot(x_sigma,volume_frac);
% xlim([0,ceil(sigma_level_max)]);
% xlabel('\sigma  level','FontSize', 15,'FontWeight','bold');
% ylabel('volume fraction','FontSize', 15,'FontWeight','bold');
% title('volume fraction over   \sigma  level','FontSize', 15,'FontWeight','bold');
% saveas(plt1,'volume_fraction_over_sigma_level.tif');
