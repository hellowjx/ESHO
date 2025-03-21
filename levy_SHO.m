function [z] = levy_SHO(pop,m,lambda)  
    num = gamma(1+lambda)*sin(pi*lambda/2); % used for Numerator 
    den = gamma((1+lambda)/2)*lambda*2^((lambda-1)/2); % used for Denominator

    sigma_u = (num/den)^(1/lambda);% Standard deviation
    u = random('Normal',0,sigma_u,pop,m);
    v = random('Normal',0,1,pop,m);
    z =u./(abs(v).^(1/lambda));
  end

