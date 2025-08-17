function [a, b, c] = my_intersect(data1, data2)
    a = intersect(data1, data2);
    b = ismember(data1, a);
    c = ismember(data2, a);
end
