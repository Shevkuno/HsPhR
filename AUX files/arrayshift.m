function [Y] = arrayshift(ARRAY)
temp1=fftshift(ARRAY,1);
Y=fftshift(temp1,2);
end

