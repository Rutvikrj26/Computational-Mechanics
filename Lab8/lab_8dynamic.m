BN = load('boundary_nodes_square1_400.txt');
Z = load('nodecoordinates_square_400.txt');
ELM = load('elem_connectivity_square1_400.txt');
nu = 0.3;
E = 1;
G = 1;
P = 0.7;
T = 7;
alpha = 0.5;
beta = 0.25;
epsilon = 0.5;
delta_t = 0.1;
l = length(Z(:,1));
u1 = zeros(l,(T/delta_t)+1);
for i = 1:l
u1(i,1) = epsilon*(Z(i,2)+0.5);
end
u2 = zeros(l,(T/delta_t)+1);
u2(:,1) = 0;
u1_dot = zeros(l,(T/delta_t)+1);
u1_dot(:,1) = 0;
u2_dot = zeros(l,(T/delta_t)+1);
u2_dot(:,1) = 0;
u1_dotdot = zeros(l,(T/delta_t)+1);
u1_dotdot(:,1) = 0;
u2_dotdot = zeros(l,(T/delta_t)+1);
u2_dotdot(:,1) = 0;
ELM_dash = ELM;
for i = 1:length(ELM(:,1))
  avg_x = (Z(ELM(i,1),2)+Z(ELM(i,2),2)+Z(ELM(i,3),2)+Z(ELM(i,4),2))/4;
  avg_y = (Z(ELM(i,1),3)+Z(ELM(i,2),3)+Z(ELM(i,3),3)+Z(ELM(i,4),3))/4;
  element = ELM(i,:);
  for j = 1:4
    if(Z(element(j),2)<avg_x && Z(element(j),3)<avg_y)
      ELM_dash(i,1) = element(j);
      end
    if(Z(element(j),2)<avg_x && Z(element(j),3)>avg_y)
      ELM_dash(i,4) = element(j);
      end
    if(Z(element(j),2)>avg_x && Z(element(j),3)<avg_y)
      ELM_dash(i,2) = element(j);
      end
    if(Z(element(j),2)>avg_x && Z(element(j),3)>avg_y)
      ELM_dash(i,3) = element(j);
      end
  endfor
endfor
Node = zeros(l,4);
for i = 1:length(ELM_dash)
  Node(ELM_dash(i,1),1) = ELM_dash(i,2);
  Node(ELM_dash(i,1),2) = ELM_dash(i,4);
  Node(ELM_dash(i,2),3) = ELM_dash(i,1);
  Node(ELM_dash(i,2),2) = ELM_dash(i,3);
  Node(ELM_dash(i,3),3) = ELM_dash(i,4);
  Node(ELM_dash(i,3),4) = ELM_dash(i,2);
  Node(ELM_dash(i,4),1) = ELM_dash(i,3);
  Node(ELM_dash(i,4),4) = ELM_dash(i,1);
end
Y = zeros(5*l,1+T/delta_t);
for i = 1:length(Node(:,1))
    if(Node(i,1)!=0&&Node(i,2)!=0&&Node(i,3)!=0&&Node(i,4)!=0)
        flag = i;
        break;
    else
        continue;
    end
end
h = (Z(Node(flag,1),2) - Z(Node(flag,3),2))/2;
b = (Z(Node(flag,2),3) - Z(Node(flag,4),3))/2;
for t = 2:(T/delta_t)+1
A = zeros(5*l,5*l);
B = zeros(5*l,1);
for i = 1:length(Node(:,1))
  if(Bound_yesno(Node(i,:))== 0)
  [x y] = Bound_check(Node(i,:));
  if(Node(i,3) == 0&&Node(i,4)==0)
      A(5*(i-1)+1,5*(i-1)+1) = -(3*E)/(2*h);
      A(5*(i-1)+1,5*(Node(i,1)-1)+1) = (4*E)/(2*h);
      A(5*(i-1)+1,5*(Node(Node(i,1),1)-1)+1) = (-E)/(2*h);
      A(5*(i-1)+1,5*(i-1)+3) = -1;
      A(5*(i-1)+1,5*(i-1)+5) = nu;
      
      A(5*(i-1)+2,5*(i-1)+2) = (-3*E)/(2*b);
      A(5*(i-1)+2,5*(Node(i,2)-1)+2) = (4*E)/(2*b);
      A(5*(i-1)+2,5*(Node(Node(i,2),2)-1)+2) = (-E)/(2*b);
      A(5*(i-1)+2,5*(i-1)+5) = -1;
      A(5*(i-1)+2,5*(i-1)+3) = nu;
      
      A(5*(i-1)+3,5*(i-1)+1) = -3*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(i,2)-1)+1) = 4*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(Node(i,2),2)-1)+1) = (-G)/(2*b);
      A(5*(i-1)+3,5*(i-1)+2) = -3*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(i,1)-1)+2) = 4*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(Node(i,1),1)-1)+2) = -(G)/(2*h);
      A(5*(i-1)+3,5*(i-1)+4) = -1;
     
   A(5*(i-1)+4,5*(i-1)+1) = 1;
   A(5*(i-1)+5,5*(i-1)+2) = 1;
  elseif(Node(i,1) == 0&&Node(i,4)==0)
      A(5*(i-1)+1,5*(i-1)+1) = (3*E)/(2*h);
      A(5*(i-1)+1,5*(Node(i,3)-1)+1) = -(4*E)/(2*h);
      A(5*(i-1)+1,5*(Node(Node(i,3),3)-1)+1) = (E)/(2*h);
      A(5*(i-1)+1,5*(i-1)+3) = -1;
      A(5*(i-1)+1,5*(i-1)+5) = nu;
      
      A(5*(i-1)+2,5*(i-1)+2) = (-3*E)/(2*b);
      A(5*(i-1)+2,5*(Node(i,2)-1)+2) = (4*E)/(2*b);
      A(5*(i-1)+2,5*(Node(Node(i,2),2)-1)+2) = (-E)/(2*b);
      A(5*(i-1)+2,5*(i-1)+5) = -1;
      A(5*(i-1)+2,5*(i-1)+3) = nu;
      
      A(5*(i-1)+3,5*(i-1)+1) = -3*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(i,2)-1)+1) = 4*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(Node(i,2),2)-1)+1) = (-G)/(2*b);
      A(5*(i-1)+3,5*(i-1)+2) = 3*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(i,3)-1)+2) = -4*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(Node(i,3),3)-1)+2) = (G)/(2*h);
      A(5*(i-1)+3,5*(i-1)+4) = -1;
      
       A(5*(i-1)+4,5*(i-1)+3) = 1;
   A(5*(i-1)+5,5*(i-1)+4) = 1;
  elseif(Node(i,1) == 0&&Node(i,2)==0)
      A(5*(i-1)+1,5*(i-1)+1) = (3*E)/(2*h);
      A(5*(i-1)+1,5*(Node(i,3)-1)+1) = -(4*E)/(2*h);
      A(5*(i-1)+1,5*(Node(Node(i,3),3)-1)+1) = (E)/(2*h);
      A(5*(i-1)+1,5*(i-1)+3) = -1;
      A(5*(i-1)+1,5*(i-1)+5) = nu;
      
      A(5*(i-1)+2,5*(i-1)+2) = (3*E)/(2*b);
      A(5*(i-1)+2,5*(Node(i,4)-1)+2) = -(4*E)/(2*b);
      A(5*(i-1)+2,5*(Node(Node(i,4),4)-1)+2) = (E)/(2*b);
      A(5*(i-1)+2,5*(i-1)+5) = -1;
      A(5*(i-1)+2,5*(i-1)+3) = nu;
      
      A(5*(i-1)+3,5*(i-1)+1) = 3*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(i,4)-1)+1) = -4*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(Node(i,4),4)-1)+1) = (G)/(2*b);
      A(5*(i-1)+3,5*(i-1)+2) = 3*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(i,3)-1)+2) = -4*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(Node(i,3),3)-1)+2) = (G)/(2*h);
      A(5*(i-1)+3,5*(i-1)+4) = -1;
      
      A(5*(i-1)+4,5*(i-1)+3) = 1;
   A(5*(i-1)+5,5*(i-1)+4) = 1;
  elseif(Node(i,3) == 0&&Node(i,2)==0)
      A(5*(i-1)+1,5*(i-1)+1) = -(3*E)/(2*h);
      A(5*(i-1)+1,5*(Node(i,1)-1)+1) = (4*E)/(2*h);
      A(5*(i-1)+1,5*(Node(Node(i,1),1)-1)+1) = (-E)/(2*h);
      A(5*(i-1)+1,5*(i-1)+3) = -1;
      A(5*(i-1)+1,5*(i-1)+5) = nu;
      
      A(5*(i-1)+2,5*(i-1)+2) = (3*E)/(2*b);
      A(5*(i-1)+2,5*(Node(i,4)-1)+2) = -(4*E)/(2*b);
      A(5*(i-1)+2,5*(Node(Node(i,4),4)-1)+2) = (E)/(2*b);
      A(5*(i-1)+2,5*(i-1)+5) = -1;
      A(5*(i-1)+2,5*(i-1)+3) = nu;
      
      A(5*(i-1)+3,5*(i-1)+1) = 3*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(i,4)-1)+1) = -4*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(Node(i,4),4)-1)+1) = (G)/(2*b);
      A(5*(i-1)+3,5*(i-1)+2) = -3*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(i,1)-1)+2) = 4*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(Node(i,1),1)-1)+2) = -(G)/(2*h);
      A(5*(i-1)+3,5*(i-1)+4) = -1;
      
      A(5*(i-1)+4,5*(i-1)+1) = 1;
   A(5*(i-1)+5,5*(i-1)+2) = 1;
  elseif(Node(i,1)==0)
      A(5*(i-1)+1,5*(i-1)+1) = (3*E)/(2*h);
      A(5*(i-1)+1,5*(Node(i,3)-1)+1) = -(4*E)/(2*h);
      A(5*(i-1)+1,5*(Node(Node(i,3),3)-1)+1) = (E)/(2*h);
      A(5*(i-1)+1,5*(i-1)+3) = -1;
      A(5*(i-1)+1,5*(i-1)+5) = nu;
      
      A(5*(i-1)+2,5*(Node(i,2)-1)+2) = E/(2*b);
        A(5*(i-1)+2,5*(Node(i,4)-1)+2) = -E/(2*b);
        A(5*(i-1)+2,5*(i-1)+5) = -1;
        A(5*(i-1)+2,5*(i-1)+3) = nu;
       
      A(5*(i-1)+3,5*(Node(i,2)-1)+1) = G/(2*b);
        A(5*(i-1)+3,5*(Node(i,4)-1)+1) = -G/(2*b); 
        A(5*(i-1)+3,5*(i-1)+2) = 3*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(i,3)-1)+2) = -4*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(Node(i,3),3)-1)+2) = (G)/(2*h);
      A(5*(i-1)+3,5*(i-1)+4) = -1;
      
      A(5*(i-1)+4,5*(i-1)+3) = 1;
   A(5*(i-1)+5,5*(i-1)+4) = 1;
  elseif(Node(i,2)==0)
       A(5*(i-1)+1,5*(Node(i,1)-1)+1) = E/(2*h);
        A(5*(i-1)+1,5*(Node(i,3)-1)+1) = -E/(2*h);
        A(5*(i-1)+1,5*(i-1)+3) = -1;
        A(5*(i-1)+1,5*(i-1)+5) = nu;
        
        A(5*(i-1)+2,5*(i-1)+2) = (3*E)/(2*b);
      A(5*(i-1)+2,5*(Node(i,4)-1)+2) = -(4*E)/(2*b);
      A(5*(i-1)+2,5*(Node(Node(i,4),4)-1)+2) = (E)/(2*b);
      A(5*(i-1)+2,5*(i-1)+5) = -1;
      A(5*(i-1)+2,5*(i-1)+3) = nu;
      
      A(5*(i-1)+3,5*(i-1)+1) = 3*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(i,4)-1)+1) = -4*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(Node(i,4),4)-1)+1) = (G)/(2*b);
      A(5*(i-1)+3,5*(Node(i,1)-1)+2) = G/(2*h);
        A(5*(i-1)+3,5*(Node(i,3)-1)+2) = -G/(2*h);
        A(5*(i-1)+3,5*(i-1)+4) = -1;    
      
    A(5*(i-1)+4,5*(i-1)+4) = 1;
   A(5*(i-1)+5,5*(i-1)+5) = 1;  
  elseif(Node(i,3)==0)
      A(5*(i-1)+1,5*(i-1)+1) = -(3*E)/(2*h);
      A(5*(i-1)+1,5*(Node(i,1)-1)+1) = (4*E)/(2*h);
      A(5*(i-1)+1,5*(Node(Node(i,1),1)-1)+1) = (-E)/(2*h);
      A(5*(i-1)+1,5*(i-1)+3) = -1;
      A(5*(i-1)+1,5*(i-1)+5) = nu;
      
      A(5*(i-1)+2,5*(Node(i,2)-1)+2) = E/(2*b);
        A(5*(i-1)+2,5*(Node(i,4)-1)+2) = -E/(2*b);
        A(5*(i-1)+2,5*(i-1)+5) = -1;
        A(5*(i-1)+2,5*(i-1)+3) = nu;
        
        A(5*(i-1)+3,5*(i-1)+2) = -3*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(i,1)-1)+2) = 4*(G)/(2*h);
      A(5*(i-1)+3,5*(Node(Node(i,1),1)-1)+2) = -(G)/(2*h);
      A(5*(i-1)+3,5*(i-1)+4) = -1;
      A(5*(i-1)+3,5*(Node(i,2)-1)+1) = G/(2*b);
        A(5*(i-1)+3,5*(Node(i,4)-1)+1) = -G/(2*b);
        
        A(5*(i-1)+4,5*(i-1)+1) = 1;
   A(5*(i-1)+5,5*(i-1)+2) = 1;
  elseif(Node(i,4)==0)
      A(5*(i-1)+1,5*(Node(i,1)-1)+1) = E/(2*h);
        A(5*(i-1)+1,5*(Node(i,3)-1)+1) = -E/(2*h);
        A(5*(i-1)+1,5*(i-1)+3) = -1;
        A(5*(i-1)+1,5*(i-1)+5) = nu;
        
        A(5*(i-1)+2,5*(i-1)+2) = (-3*E)/(2*b);
      A(5*(i-1)+2,5*(Node(i,2)-1)+2) = (4*E)/(2*b);
      A(5*(i-1)+2,5*(Node(Node(i,2),2)-1)+2) = (-E)/(2*b);
      A(5*(i-1)+2,5*(i-1)+5) = -1;
      A(5*(i-1)+2,5*(i-1)+3) = nu;
      
      A(5*(i-1)+3,5*(Node(i,1)-1)+2) = G/(2*h);
        A(5*(i-1)+3,5*(Node(i,3)-1)+2) = -G/(2*h);
        A(5*(i-1)+3,5*(i-1)+1) = -3*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(i,2)-1)+1) = 4*(G)/(2*b);
      A(5*(i-1)+3,5*(Node(Node(i,2),2)-1)+1) = (-G)/(2*b);
      A(5*(i-1)+3,5*(i-1)+4) = -1;
      
      A(5*(i-1)+4,5*(i-1)+4) = 1;
   A(5*(i-1)+5,5*(i-1)+5) = 1;
 end  
  else
        A(5*(i-1)+1,5*(Node(i,1)-1)+1) = E/(2*h);
        A(5*(i-1)+1,5*(Node(i,3)-1)+1) = -E/(2*h);
        A(5*(i-1)+1,5*(i-1)+3) = -1;
        A(5*(i-1)+1,5*(i-1)+5) = nu;
      
        A(5*(i-1)+2,5*(Node(i,2)-1)+2) = E/(2*b);
        A(5*(i-1)+2,5*(Node(i,4)-1)+2) = -E/(2*b);
        A(5*(i-1)+2,5*(i-1)+5) = -1;
        A(5*(i-1)+2,5*(i-1)+3) = nu;
      
        A(5*(i-1)+3,5*(Node(i,2)-1)+1) = G/(2*b);
        A(5*(i-1)+3,5*(Node(i,4)-1)+1) = -G/(2*b);
        A(5*(i-1)+3,5*(Node(i,1)-1)+2) = G/(2*h);
        A(5*(i-1)+3,5*(Node(i,3)-1)+2) = -G/(2*h);
        A(5*(i-1)+3,5*(i-1)+4) = -1;
        
        A(5*(i-1)+4,5*(Node(i,1)-1)+3) = 1/(2*h);
        A(5*(i-1)+4,5*(Node(i,3)-1)+3) = -1/(2*h);
        A(5*(i-1)+4,5*(Node(i,2)-1)+4) = 1/(2*b);
        A(5*(i-1)+4,5*(Node(i,4)-1)+4) = -1/(2*b);
        A(5*(i-1)+4,5*(i-1)+1) = -1/(beta*delta_t^2);
        B(5*(i-1)+4) = u1(i,t-1)/(beta*delta_t^2)+u1_dot(i,t-1)/(beta*delta_t)+(1-2*beta)*u1_dotdot(i,t-1)/(2*beta);
        
        
        A(5*(i-1)+5,5*(Node(i,1)-1)+4) = 1/(2*h);
        A(5*(i-1)+5,5*(Node(i,3)-1)+4) = -1/(2*h);
        A(5*(i-1)+5,5*(Node(i,2)-1)+5) = 1/(2*b);
        A(5*(i-1)+5,5*(Node(i,4)-1)+5) = -1/(2*b);
        A(5*(i-1)+5,5*(i-1)+2) = -1/(beta*delta_t^2);
        B(5*(i-1)+5) = u2(i,t-1)/(beta*delta_t^2)+u2_dot(i,t-1)/(beta*delta_t)+(1-2*beta)*u2_dotdot(i,t-1)/(2*beta);
   end
end
if(t==2)
Bdash = B;
end
Y(:,t) = -A\B;
   for i = 1:l
   u1(i,t) = Y(5*(i-1)+1,t);
   u2(i,t) = Y(5*(i-1)+2,t);
   u1_dotdot(i,t) = (1/(2*beta))*((2/delta_t^2)*(u1(i,t)-u1(i,t-1)-delta_t*u1_dot(i,t-1))-(1-2*beta)*(u1_dotdot(i,t-1)));
   u1_dot(i,t) = u1_dot(i,t-1)+delta_t*((1-alpha)*u1_dotdot(i,t-1)+alpha*u1_dotdot(i,t));
   u2_dotdot(i,t) = (1/(2*beta))*((2/delta_t^2)*(u2(i,t)-u2(i,t-1)-delta_t*u2_dot(i,t-1))-(1-2*beta)*(u2_dotdot(i,t-1)));
   u2_dot(i,t) = u2_dot(i,t-1)+delta_t*((1-alpha)*u2_dotdot(i,t-1)+alpha*u2_dotdot(i,t));
   end
endfor
#while(true)
for j = 1:(T/delta_t)+1
  final = zeros(l,3);
for i = 1:l
  final(i,2) = Z(i,2) + u1(i,j);
  final(i,3) = Z(i,3) + u2(i,j);
endfor
plot(final(:,2),final(:,3),'*');
axis([-1 1.5 -1 1]);
if(j==1)
pause(1);
else
pause(0.001)
end
end
#pause(1);
#end