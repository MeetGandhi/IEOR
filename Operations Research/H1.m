clear all
Z=0;
IT=0;
Z1=0;
c=input('Input Cost Coefficient Array: ');
c=transpose(c);
l_c=length(c);
x=zeros(l_c,1);
b=input('Input Constraint Constants as an Array: ');
b=transpose(b);
B=b;
a=input('Input Constraint Coefficients as a Matrix: ');
l_b=length(b);
x_s=zeros(l_b,1);
z=zeros(l_b,1);
I=eye(l_b);
L=[1,transpose(-1*c),transpose(z);z,a,I];
R=[Z;x;x_s];
Q=[Z1;b];
E=[L,Q];
while min(L(1,:))~=0
    IT=IT+1;
    disp(strcat('Iteration: ',num2str(IT)))
    [m,n]=find(L==min(L(1,:)));
    B_t=[];
    B_tt=[];
    for i=1:l_b
        b_t=b(i)/L(i+1,n);
        B_tt=[B_tt,b_t];
        if L(i+1,n)~=0
            b_t=b(i)/L(i+1,n);
            if b_t>0
                B_t=[B_t,b_t];
            end
        end
    end
    [o,p]=find(B_tt==min(B_t));
    if length(o)~=1 || length(p)~=1
        o=o(1);
        p=p(1);
    end
    p=p+1;
    [l_r,l_c]=size(E);
    for j=1:l_c
        E(p,j)=E(p,j)/L(p,n);
    end
    
    for i=1:l_r
        if i~=p
            mi=-1*L(i,n);
            for j=1:l_c
                E(i,j)=E(i,j)+(mi*E(p,j));
            end
        end
    end
    bv=[];
    nbv=[];
    bfs=[];
    for i=2:l_c
        if norm(E(:,i))==1.0
            bv=[bv,i-1];
            for j=2:l_r
                if E(j,i)==1 || E(j,i)==1.0 
                    bfs=[bfs,E(j,l_c)];
                end
            end
        else
            nbv=[nbv,i-1];
        end
    end
   
    disp(strcat('Basic Variables: ',num2str(bv)))
    disp(strcat('Current BFS point (Values of BVs): ',num2str(bfs)))
    disp(strcat('Non Basic Variables: ',num2str(nbv)))
    disp(strcat('Z=',num2str(E(1,l_c))))
    L=E;
    b=E(2:l_r,l_c);
end
disp(strcat('Shadow Prices of consecutive constraints in order: ', num2str(E(1,length(x)+2:l_c-1))))