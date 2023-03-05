
global SUCCESS FAIL
SUCCESS = 0; 
FAIL = 1;

diary gaussElimOutput

A = [3,1,2;6,3,4;3,1,5];
b = [0;1;3];
disp("A = ")
disp(A)
disp("b = ")
disp(b)
fprintf("\r\n")
luSolve(A,b);
fprintf("\n\n")


C = [4,2,0;4,4,2;2,2,3];
d = [2;4;6];
disp("A = ")
disp(C)
disp("b = ")
disp(d)
fprintf("\n")
luSolve(C,d);
fprintf("\n\n")

diary off

function[L,U,s]=luFactor(a)
    global SUCCESS FAIL
    s=SUCCESS;
    [~,n]=size(a);
    
    L=zeros(n,n);
    
    U=zeros(n,n);
    
    for i=1:n
        if abs(a(i,i))<0
            error('zero pivot encountered');
            s = FAIL;
            break
        end

        for k=1:i-1

            L(i,k)=a(i,k);
            for j=1:k-1
                L(i,k)=L(i,k)-L(i,j)*U(j,k);
            end

            L(i,k)=L(i,k)/U(k,k);
        end

        for k=1:n
            U(i,k)=a(i,k);
            for j=1:i-1
                U(i,k)=U(i,k)-L(i,j)*U(j,k);
            end
        end
        L(i,i)=1;

    end
    
end

function [x]=luFactorSolve(L,U,b)
    [n,~]=size(b);
    c=zeros(n,1);
    c(1)=b(1)/L(1,1);
    for k = 2:n
        s=0;
        for j = 1:k-1
            s = s + L(k,j)*c(j);
        end
        c(k) = b(k) - s;
    end
    x = zeros(n,1);
    for k = n:-1:1
        t=0;
        for j = k+1:n
            t = t + U(k,j)*x(j);
        end
        x(k) = (c(k) - t)/U(k,k);
    end

end

function [x]=luSolve(a,b)
    global SUCCESS FAIL
    [L,U,s] = luFactor(a);
    if s == SUCCESS
        disp('L =')
        disp(L)
        disp('U =')
        disp(U)
    
        x = luFactorSolve(L,U,b);
        disp('x =')
        disp(x)
    end
end









