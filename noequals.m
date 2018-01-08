function [dataout] = noequals(data)
x = data;
while length(data)~=length(unique(x))
    for  i =2:1:length(x)
        if x(i) <= x(i-1)
            x(i) = x(i)+0.0001;
        end
    end
end
dataout = x;
end