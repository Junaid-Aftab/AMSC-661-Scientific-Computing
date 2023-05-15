function FreeSchrodinger()
close all
fsz = 15; % fontsize

% solves u_t = i/2 u_{xx};

%Space Parameters
dx = 0.01;
L = 40;                 %Length of x-interval
N = round(L/dx);        %Coallocation Points
x = linspace(-L/2,L/2,N+1);
x(end) = [];

% Time Parameters
dt = 0.01;
tmax = 0.4;
t = 0;
pause_duration = 0.05;
epsilon = 1e-15;

% Frequencies
k = -N/2 : (N/2 - 1); % wave numbers
freq = k.*(2*pi/L); 
freq2 = freq.^2;

% Initial Data
k0 = 10;
sigma = 0.1;
u0 = power( 1/(2*pi*sigma^2), 1/4).*exp(- x.^2/(4*sigma^2) + 1i*k0.*x);

% Solution via two methods
u_dft = u0;
u_rk = u0;

% Figure
figure(1); clf;   % Create a new figure
hold on;
legendLabels = cell(1, 6);  % Initialize cell array to store legend labels

while(t <= tmax+epsilon)

    %METHOD 1: Discrete Fourier Transform

    % u0 in Fourier space
    v_dft = fftshift(fft(u_dft));

    % solve in Fourier space.
    vnew_dft = exp(-0.5i*dt.*freq2).*v_dft;

    % return to real space
    unew_dft = ifft(ifftshift(vnew_dft)); 

    %METHOD 2: Method of Lines + RK4

    % RK4 step
    k1=rhs(dx,u_rk);
    k2=rhs(dx,u_rk+0.5*dt*k1);
    k3=rhs(dx,u_rk+0.5*dt*k2);
    k4=rhs(dx,u_rk+dt*k3);
    unew_rk=u_rk+(dt/6).*(k1+2*k2+2*k3+k4);

    %update plots
    if abs(5*t/tmax-round(5*t/tmax)) < 10*epsilon
        plot(x, abs(u_dft),'LineWidth',2);
        hold on;
        legendLabels{int32(round(5*t/tmax))+1} = sprintf('Time = %g', t);
        xlabel('x','FontSize',fsz);
        ylabel('Norm Squared of Wave Function','FontSize',fsz);
        set(gca,'FontSize',fsz);
        grid on;
        drawnow
    end
    %disp(legendLabels)

    % Update varibales etc.
    probability = prob(sigma,k0,x,t);
    trapz(x,probability);
    u_dft = unew_dft;
    u_rk = unew_rk;
    t=t+dt;

    %set(hpic1,'xdata',x,'ydata', sqrt(probability) );
    %set(hpic2,'xdata',x,'ydata', abs(u_dft) );
    
    pause(pause_duration);
end

% Legend
legend(legendLabels, 'Location', 'best');

function rhs = prob(sigma,k,position,time)
    A1 = power(2/pi, 1/4);
    A2 = power( (1i*time)/sigma + 2*sigma, 1/2);
    A3 = 1i*position.^2 -2*k*sigma^2*(k*time - 2.*position);
    A4 = 2*(time - 2i*sigma^2);
    wave = (A1/A2)*exp(A3/A4);
    wave_conj = conj(wave);
    rhs=wave.*wave_conj;

function RHS=rhs(dx,u_inter)
    shiftedLeft = circshift(u_inter, [0, -1]);
    shiftedRight = circshift(u_inter, [0, +1]);
    RHS = (1i/(2*dx^2)).*( shiftedLeft + shiftedRight - 2.*u_inter );


%Plot initial figures
%figure(1); clf; 
%probability = prob(sigma,k0,x,t);
%hpic1 = plot(x, sqrt(probability),'LineWidth',2,'color','b'); 
%hold on;
%hpic2 = plot(x, abs(u0),'LineWidth',2,'color','r'); 
%xlabel('x','FontSize',fsz);
%ylabel('Norm Squared of Wave Function','FontSize',fsz);
%set(gca,'FontSize',fsz);
%grid
%drawnow
