function [Y] = arrayIshift(ARRAY)
Y=ifftshift(ifftshift(ARRAY,1),2);
end

