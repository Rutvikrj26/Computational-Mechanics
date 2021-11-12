n = input("value of n: ");
m = input("Value of m: ");
A=zeros(n,n);
for i=1:m
    for j= 1:n+1-i
        A(j,i+j-1)=rand();
        if(i~=1)
            A(i+j-1,j)=rand();
        end
    end
end
B=A;
L=eye(n,n);
for b=1:n-1
     g=min(n,b+m-1);
     for k=b+1:g
        L(k,b)=A(k,b)/A(b,b);
     end
     for h=b+1:g
     A(h,[1:g])=A(h,[1:g])-A(b,[1:g])*L(h,b);
     end
end
disp(L);
disp('');
disp(A);
disp('');
disp(B);
disp('');
disp(L*A);