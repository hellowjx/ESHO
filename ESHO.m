function [BestFitness,BestPosition,Convergence_curve]=ESHO(N,Max_iter,lb,ub,dim,fobj)


%Logistic Chaos Mapping Initialzation
x1(1)=rand;
for i=2:N*dim
    x1(i)=4.*x1(i-1).*(1-x1(i-1));
end

k=1;
if size(ub,2)>1
    Sea_horses=zeros(N,dim);
    for i = 1:N
        for j=1:dim
            Sea_horses(i,j)=x1(k).*(ub(j)-lb(j))+lb(j);
            k = k+1;
        end
    end
else
    Sea_horses=zeros(N,dim);
    for i = 1:N
        for j=1:dim
            Sea_horses(i,j)=x1(k).*(ub-lb)+lb;
            k = k+1;
        end
    end
end



Sea_horsesFitness = zeros(1,N);  %Fitness
Convergence_curve=zeros(1,Max_iter); %Convergence_curve


%%Calculate fitness
for i=1:N
    Sea_horsesFitness(1,i)=fobj(Sea_horses(i,:));
end

%%find the best individual
[~,sorted_indexes]=sort(Sea_horsesFitness);
BestPosition=Sea_horses(sorted_indexes(1),:);
BestFitness = Sea_horsesFitness(sorted_indexes(1));
Convergence_curve(1)=BestFitness;

t=1;
u=0.05;  
v=0.05;
l=0.05;
trial=zeros(N,dim);


while t<Max_iter
    beta=randn(N,dim);  
    Elite=repmat(BestPosition,N,1);
    r1=randn(1,N);
    Step_length=levy_SHO(N,dim,1.5);
    %Mobile
    for i=1:N  
        for j=1:dim
            if r1(i)>0
                r=rand();
                theta=r*2*pi;
                row=u*exp(theta*v);
                x=row*cos(theta);
                y=row*sin(theta);
                z=row*theta;
                Sea_horses_new1(i,j)=Sea_horses(i,j)+Step_length(i,j)*((Elite(i,j)-Sea_horses(i,j)).*x.*y.*z+Elite(i,j));%Eq.(2)
            else
                Sea_horses_new1(i,j)=Sea_horses(i,j)+rand()*l*beta(i,j)*(Sea_horses(i,j)-beta(i,j)* Elite(i,j));%Eq.(2)
            end
        end
    end

    for i=1:N
        %Bound the variable
        Tp=Sea_horses_new1(i,:)>ub;
        Tm=Sea_horses_new1(i,:)<lb;
        Sea_horses_new1(i,:)=(Sea_horses_new1(i,:).*(~(Tp+Tm)))+ub.*Tp+lb.*Tm;
    end

  
    %Predatory
    for i=1:N
        for j=1:dim
            r2(i)=rand();
            alpha=(1-t/Max_iter)^(2*t/Max_iter);
            if r2(i)>=0.1
                Sea_horses_new2(i,j)=alpha*(Elite(i,j)-rand()*Sea_horses_new1(i,j))+(1-alpha)*Elite(i,j);  %Eq.(6)
            else
                Sea_horses_new2(i,j)=(1-alpha)*(Sea_horses_new1(i,j)-rand()*Elite(i,j))+alpha*Sea_horses_new1(i,j);  %Eq.(6)
            end
        end
    end

    for i=1:N 
        %Bound the variable
        Tp=Sea_horses_new2(i,:)>ub;
        Tm=Sea_horses_new2(i,:)<lb;
        Sea_horses_new2(i,:)=(Sea_horses_new2(i,:).*(~(Tp+Tm)))+ub.*Tp+lb.*Tm;
        %Calculated all sea horses' fitness values
        Sea_horsesFitness1(i)=fobj(Sea_horses_new2(i,:));
    end



    c1 = 0.2+t*((0.5)/Max_iter);
    c2 = 0.2+t*((0.5)/Max_iter);
    c3 = 0.7-t*((0.5)/Max_iter);
    c4 = 0.7-t*((0.5)/Max_iter);
    f1=0.1;f2=0.1;f3=0.1;
    restart = 1+t*((3)/Max_iter);
    
    %Adaptive Multiple mutation strategy(AMMS)
    for i = 1:N
        X = Sea_horses_new2;
        [~,sortIndex]=sort(Sea_horsesFitness1);
        best_index = sortIndex(1);

    
        v1=X(i,:); v2=X(i,:); v3=X(i,:); v4=X(i,:); 

        if rand() < c1   
            for j=1:dim
                r1=randi(N);r2=randi(N);r3=randi(N);
                v1(1,j) = X(r1,j) + f1*(X(best_index,j)-X(r1,j)) + f1*(X(r2,j)-X(r3,j));             % DE/rand-to-best/1
            end
            Flag4ub=v1(1,:)>ub;
            Flag4lb=v1(1,:)<lb;
            v1(1,:)=(v1(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        end

        if rand() < c2
            for j=1:dim
                r4=randi(N);r5=randi(N);r6=randi(N);r7=randi(N);
                v2(1,j) = X(best_index,j) +f2*(X(r4,j)-X(r5,j))+ f2*(X(r6,j)-X(r7,j));              % DE/best/2
            end
            Flag4ub=v2(1,:)>ub;
            Flag4lb=v2(1,:)<lb;
            v2(1,:)=(v2(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        end


        if rand() < c3
            for j=1:dim
                r8=randi(N);r9=randi(N);r10=randi(N);
                v3(1, j) =  X(i,j)  + f3*(X(r8,j)-X(i,j))+f3*(X(r9,j)-X(r10,j));                   % DE/current to rand/1
            end
            Flag4ub=v3(1,:)>ub;
            Flag4lb=v3(1,:)<lb;
            v3(1,:)=(v3(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        end



        if rand() < c4
            for j=1:dim
                r11=randi(N);
                v4(1,j) = X(i,j)+CauchyRand(0,1)*X(r11,j);
            end

            Flag4ub=v4(1,:)>ub;
            Flag4lb=v4(1,:)<lb;
            v4(1,:)=(v4(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        end



        fv1 = fobj(v1);   fv2 = fobj(v2);   fv3 = fobj(v3);    fv4 = fobj(v4);
        newPosition = [v1;v2;v3;v4];
        newFitness = [fv1,fv2,fv3,fv4];
        [~,newIndex]=sort(newFitness);

        if newFitness(newIndex(1)) < Sea_horsesFitness1(i)
            X(i,:) = newPosition(newIndex(1),:);
            Sea_horses_new2(i,:) = X(i,:);
            Sea_horsesFitness1(i)=newFitness(newIndex(1)) ;
            trial(i) = 0;
        else
            trial(i) = trial(i)+ 1;
        end
        
        %Adaptive Restart Mechanism(ARM)
        if trial(i) > restart
            t1 = zeros(1,dim); t2 = zeros(1,dim);
            t1(1,:) = (ub-lb)*rand+lb;
            t2(1,:) = (ub+lb)*rand-X(i,:);
            Flag4ub=t2(1,:)>ub;
            Flag4lb=t2(1,:)<lb;
            t2(1,:)=(t2(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            if fobj(t1) < fobj(t2)
                X(i,:) = t1(1,:);
                Sea_horses_new2(i,:) = X(i,:);
                Sea_horsesFitness1(i)=fobj(t1) ;
            else
                X(i,:) = t2(1,:);
                Sea_horses_new2(i,:) = X(i,:);
                Sea_horsesFitness1(i)=fobj(t2) ;
            end
            trial(i) = 0;
        end
    end %end for



    [~,sorted_indexes]=sort(Sea_horsesFitness1);
    %Select the best pop and sort it
    Sea_horses=Sea_horses_new2(sorted_indexes(1:N),:);
    SortfitbestN = Sea_horsesFitness1(sorted_indexes(1:N));

    %Update the optimal solution
    if SortfitbestN(1)<BestFitness
        BestPosition=Sea_horses(1,:);
        BestFitness=SortfitbestN(1);
    end



    Convergence_curve(t)=BestFitness;
    t = t + 1;


end
end

function [z] = levy_SHO(pop,m,lambda)  
    num = gamma(1+lambda)*sin(pi*lambda/2); % used for Numerator 
    den = gamma((1+lambda)/2)*lambda*2^((lambda-1)/2); % used for Denominator

    sigma_u = (num/den)^(1/lambda);% Standard deviation
    u = random('Normal',0,sigma_u,pop,m);
    v = random('Normal',0,1,pop,m);
    z =u./(abs(v).^(1/lambda));
end

function [cauchy] = CauchyRand(m,c)
cauchy = c*tan(pi*(rand()-0.5)) + m;
end


