clc; clear; clf;

% SET MODEL

c=1;                               %wave velocity
Nx=100;                            %numper of points
Nt=60;                             %whole time
Lx=1;                             %length in 1D 
dx=Lx/(Nx-1);                      %step on x axes  
dt=[(0.9*dx)/c,(dx)/c,(1.1*dx)/c]; %time step -> 0.9dx/c
                                   %          -> 1.0dx/c (magic step        
                                   %          -> 1.1dx/c
                                   
% EXPLICIT=========================



for Dt=1:3
     
    
r=c*(dt(Dt)/dx);
    u_1=zeros(Nx,1);                %t+dt  
    u_o=zeros(Nx,1);                %t
    u_m1=zeros(Nx,1);               %t-dt       
    
    u_time=zeros(Nx,5);             %use array to hold on u(n+1) every 10  
                                    %time steps
    u_m1(2:11,:)=1;                 %initial condition for t-dt
    u_o(3:12,:)=1;                  %initial contition for t using step 
                                    %to the right 
   time=1;
    for nt=1:Nt
         for j=2:Nx-1               % calculating EXPLICIT t+dt  array
           u_1(j,1)=r^2.*(u_o(j+1,1)-2.*u_o(j,1)+u_o(j-1,1))+2.*u_o(j,1)-u_m1(j,1);
         end
       
            
        if nt==10*time
           u_time(:,time)= u_1;   %hold on the solution every 10 time steps
            time=time+1;          %EXPLICIT
         end
    
        u_m1=u_o;
        u_o=u_1;   
    end 

    %IMPLICIT=========================
    
    bo=0.25;
    uI_1=zeros(Nx,1);         
    uI_o=zeros(Nx,1);         
    uI_m1=zeros(Nx,1);                      
    
  
    
    

% define tridiagonial array L    
a = 2;
b = -1;
c = -1;
L = diag(a*ones(1,Nx)) + diag(b*ones(1,Nx-1),1) + diag(c*ones(1,Nx-1),-1);

    I=eye(Nx);                  %define indentity array 
    uI_time=zeros(Nx,5);
    A=bo*L+r^2*I;
    
    B=((2*bo-1)/2)*L+r^2*I;
    
    uI_m1(2:11,:)=1;
    uI_o(3:12,:)=1;
    timeI=1;
    for j=1:Nt
        uI_1= 2*(A\B) * uI_o - uI_m1 ; %calculating solution array IMPLICIT    
        if j==10*timeI
            uI_time(:,timeI)= uI_1; 
            timeI=timeI+1;
        end
        uI_m1=uI_o;
        uI_o=uI_1; 
    end 
    
st=[0.9,1.0,1.1];
figure(Dt)                       %figures EXPLICIT-IMPLICIT for every dt
a=st(Dt);
sgt = sgtitle(['Explicit  dt=',num2str(a),'dx/c Implicit'],'Color','red');
sgt.FontSize = 20;
 

    subplot(5,2,1)
    x=linspace(0,Lx,Nx);
    col = 'b'; LW = 1.2;
    plot(x,u_time(:,2),col,'linewidth',LW)
    grid on; box on;
    title('t=20 dt')
    
    subplot(5,2,2)
    x=linspace(0,Lx,Nx);
    col = 'b'; LW = 1.2;
    plot(x,uI_time(:,2),col,'linewidth',LW)
    grid on; box on;
    title('t=20')
    
    
    subplot(5,2,3)
    x=linspace(0,Lx,Nx);
    plot(x,u_time(:,3),col,'linewidth',LW)
    grid on; box on;
    title('t=30')
    
    subplot(5,2,4)
    x=linspace(0,Lx,Nx);
    plot(x,uI_time(:,3),col,'linewidth',LW)
    grid on; box on;
    title('t=30')
    
    
    
    subplot(5,2,5)
    x=linspace(0,Lx,Nx);
    plot(x,u_time(:,4),col,'linewidth',LW)
    grid on; box on;
    title('t=40')
    ylabel('U(x,t)')
    
    subplot(5,2,6)
    x=linspace(0,Lx,Nx);
    plot(x,uI_time(:,4),col,'linewidth',LW)
    grid on; box on;
    title('t=40')
    ylabel('U(x,t)')
    
    subplot(5,2,7)
    x=linspace(0,Lx,Nx);
    plot(x,u_time(:,5),col,'linewidth',LW)
    grid on; box on;
    title('t=50')
    
    subplot(5,2,8)
    x=linspace(0,Lx,Nx);
    plot(x,uI_time(:,5),col,'linewidth',LW)
    grid on; box on;
    title('t=50')
    
    
    subplot(5,2,9)
    x=linspace(0,Lx,Nx);
    plot(x,u_time(:,6),col,'linewidth',LW)
    grid on; box on;
    title('t=60')
    xlabel('x');
    
    subplot(5,2,10)
    x=linspace(0,Lx,Nx);
    plot(x,uI_time(:,6),col,'linewidth',LW)
    grid on; box on;
    title('t=60')
    xlabel('x');
    
   
end