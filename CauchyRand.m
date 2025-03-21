function [cauchy] = CauchyRand(m,c)
cauchy = c*tan(pi*(rand()-0.5)) + m;
end

