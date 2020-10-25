RGB = imread('BadQuality1.png');
I=rgb2gray(RGB);
img = I(11:20,2:11);
x = 0:7; x= x';
result = graycomatrix(img,'Offset',[1 0]);
p_x_y = result/sum(result(:));
p_x = sum(p_x_y, 1); p_x = p_x';
p_y = sum(p_x_y, 2);
mu_x = sum(x.* p_x);
mu_x_squared = sum(x.^2 .* p_x);
var_x = mu_x_squared - mu_x^2;

mu_x_times_y = sum(sum((x*x').*p_x_y));
correlation = (mu_x_times_y - mu_x*mu_x)/var_x;
stats = graycoprops(result,{'Correlation'})