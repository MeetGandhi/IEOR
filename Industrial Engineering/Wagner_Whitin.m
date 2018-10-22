function [ xsolution,totalcost ] = Wagner_Whitin(demand,n,sc,h )
% Solve the time varying demand problem with Wagner-Whitin algorithm 
% Written by W.Boonphakdee,D.Eng., Oct 3,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of input parameter
% demand=row vector of demand
% n=number of period
% sc=setup cost
% h=holding cost
% Description of output parameter
% xsolution= row vector of replenishment (Lot size)
% totalcost= total cost of setup and holding cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If users have any questions you can contact me : warutboon@yahoo.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% totaldemand=sum(demand);
%[m,n]=size(demand);% n = number of period
% sc=54;% setup cost
% h=0.4;% holding cost

%% The functional equation matrix and it functional equation of the last period 
stage=1;
% Construct the functional equation matrix
fnmat1=zeros(stage+1,stage+1);% r1(S1,d1) matrix of stage 1
xmat=zeros(stage+1,stage);% dn matrix
% Determine the number of inventory and the number of replenishment
stock=zeros(stage+1,stage);
replen=zeros(stage,stage+1);
for i=1:stage+1
    if i==1
       stock(i,stage)=0;
       replen(i,stage)=0;
    else
       stock(i,stage)=demand(n-(stage-1));
       replen(stage,i)=stock(i,stage);
    end
end
% The functional equation matrix
for i=1:stage;
    for j=1:stage+1;
        if  replen(i,j)>0
            fnmat1(i,j)= sc+h*(stock(i)+replen(i,j)-demand(n));
        end
    end
end
 % Minimal the functional equation
 minfn=zeros(n+1,stage);
 minfnadd=zeros(stage+1,stage);
 R=inf;
 for i=1:stage+1
     for j=1:stage+1
         if fnmat1(i,j)>0 & fnmat1(i,j)<=R
             R=fnmat1(i,j);
             minfn(i,stage)=fnmat1(i,j);
             xmat(i,stage)=replen(i,j);
             minfnadd(i,stage)=fnmat1(i,j);
         end
     end
 end
 %% The functional equation matrix and it functional equation of the other period 
 
 for stage=2:n
   fnmat2=zeros(stage+1,stage+1);  
   if stage <n
    % Determine the number of inventory and the number of replenishment
     countstock=0;
     sumstock=0;
     for i=1:stage+1
        if i==1
           stock(i,stage)=0;
           replen(stage,i)=0;
        else
           countstock=countstock+1;
           sumstock=sumstock+demand(n-stage+countstock);
           stock(i,stage)=sumstock;
           replen(stage,i)=sumstock;
        end
     end 
     % The functional equation matrix
       countreplen1=0;
       countreplen2=0;
       for i=1:stage+1;
          for j=1:stage+1;
              if  i==1 & replen(stage,j)>0
                  countreplen1=countreplen1+1;
                  fnmat2(i,j)= sc+h*(stock(i,stage)+replen(stage,j)-demand(n-stage+1))+minfn(countreplen1,stage-1);
                  fnadd(i,j)=minfn(countreplen1,stage-1);
              else
                  if  stock(i,stage)>0 & replen(stage,j)==0
                      countreplen2=countreplen2+1;
                      fnmat2(i,j)=h*(stock(i,stage)-demand(n-stage+1))+minfn(countreplen2,stage-1);
                      fnadd(i,j)=minfn(countreplen2,stage-1);
                  end
              end
          end
       end 
    % Searching the minimal functional equation
     
      for i=1:stage+1
           R=inf;
          for j=1:stage+1
              if fnmat2(i,j)>0 & fnmat2(i,j)<=R
                 R=fnmat2(i,j);
                 minfn(i,stage)=fnmat2(i,j);
                 xmat(i,stage)=replen(stage,j);
              end
          end
      end
   else % stage =n
      % Transform the previous fnmat matrix to null matrix
      for i=1:stage+1
          for j=1:stage+1
              fnadd(i,j)=0;
          end
      end
     % Determine the number of inventory and the number of replenishment
     ncountstock=0;
     nsumstock=0;
     for i=1:1
         for j=1:n
        
           ncountstock=ncountstock+1;
           nsumstock=nsumstock+demand(j);
           stock(j,stage)=0;
           replen(stage,j+1)=nsumstock;
         end
     end 
     % The functional equation matrix
       ncountreplen1=0;
       ncountreplen2=0;
       for i=1:stage+1
          for j=1:stage+1;
              if i==1 & replen(stage,j)>0
                  
                  ncountreplen1=ncountreplen1+1;
                  fnmat2(i,j)= sc+h*(stock(i,stage)+replen(stage,j)-demand(n-stage+1))+minfn(ncountreplen1,stage-1);
                  fnadd(i,j)=minfn(ncountreplen1,stage-1);
              
              end
          end
       end 
    % Searching the minimal functional equation
      R=inf;
      for i=1:stage+1
          for j=1:stage+1
              if fnmat2(i,j)>0 & fnmat2(i,j)<=R
                 R=fnmat2(i,j);
                 minfn(i,stage)=fnmat2(i,j);
                 xmat(i,stage)=replen(stage,j);
              end
          end
      end
   end
 end
  %% Optimal solution
  xsolution=zeros(1,stage);
  remaindemand=zeros(1,stage);
  countstage=1;
  remaindemand(1)=xmat(1,n)-demand(1);
  xsolution(1)=xmat(1,n);
  for stage=n-1:-1:1
    
      if demand(n-stage+1)<remaindemand(countstage)
          countstage=countstage+1;
          xsolution(countstage)=0;
          remaindemand(countstage)=remaindemand(n-stage)-demand(n-stage+1);
      elseif demand(n-stage+1)==remaindemand(n-stage)
          countstage=countstage+1;
          xsolution(countstage)=0;
          remaindemand(n-stage+1)=0;
      elseif demand(n-stage+1)==xmat(1,n-countstage)
          countstage=countstage+1;
          remaindemand(n-stage+1)=0;
          xsolution(countstage)=xmat(1,n-countstage+1);
      elseif demand(n-stage+1)<xmat(1,n-countstage)
          
          remaindemand(n-stage+1)=xmat(1,n-countstage)-demand(n-stage+1);
         countstage=countstage+1; 
          xsolution(countstage)=xmat(1,n-countstage+1);
      end
  
  end
totalcost=minfn(1,n);
 


end

