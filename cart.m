%Cartesian products of two arrays X and Y
function xy = cart(x, y)
xy=[];
for i = x
  for j = y 
    xy = [xy, [i;j] ]; 
    end
end
end