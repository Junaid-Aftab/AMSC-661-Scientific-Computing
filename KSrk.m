function KSrk()
close all
fsz = 20; % fontsize
% solves u_t = - u_{xx} - u_{xxxx} - (0.5u^2)_x

N = 256;
L = 32*pi;
x = linspace(-L/2,L/2,N+1);
x(end) = [];
k = -N/2 : (N/2 - 1); % wave numbers

% initial data
u0=cos(x/16).*(1+sin(x/16));

dt = 0.1; % time step

figure(1); clf; 
hpic = plot(x,u0,'LineWidth',2,'color','r'); % plot the initial condition for the solution to be computed
hold on;
xlabel('x','FontSize',fsz);
ylabel('u','FontSize',fsz);
set(gca,'FontSize',fsz);
grid
drawnow
pad=0.1*(max(abs(u0)));
xlim([-L/2 L/2 ]);

tmax = 200;
t = 0;
u_all = zeros(N,tmax/dt+1);

freq = k.*(2*pi/L); % frequencies
freq2 = freq.^2;
freq4 = freq2.^2;
e3=exp((freq2 - freq4)*dt); % in the Fourier space, uhat = e3.*vhat

while (t<tmax) 
    t=t+dt;
    vhat=fftshift(fft(u0)); % v in the Fourier space
    % RK4 step in the Fourier space
    k1=rhs(0,vhat);
    k2=rhs(0.5*dt,vhat+0.5*dt*k1);
    k3=rhs(0.5*dt,vhat+0.5*dt*k2);
    k4=rhs(dt,vhat+dt*k3);
    vhat_new=vhat+dt*(k1+2*k2+2*k3+k4)/6;
    % return to the original space and the original variable u
    unew=ifft(ifftshift(e3.*vhat_new)); % return to u in the x-space
    u0=unew;
    set(hpic,'xdata',x,'ydata',real(unew));
    %drawnow
    u_all(:,floor(t/dt)) = transpose(real(u0));
end

% Plot the surface
figure;
[T, X] = meshgrid( 0:dt:tmax, -L/2:L/(N-1) :L/2 );
surf(X, T, u_all, 'EdgeColor', 'none');
title('Surface Plot');
xlim([0, x(end)])
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
 
% Apply color mapping using imagesc
colormap('jet');  % Choose a color map
shading interp;  % Interpolate colors between vertices
colorbar;  % Display color bar
 
% Step 2: Get coordinates and color data
coords = [T(:), X(:), u_all(:)];  % Combine X, Y, and Z coordinates
colors = u_all(:);               % Use Z values as color data
 
% Step 3: Project coordinates onto a 2D plane
proj_coords = coords(:, 1:2); % Use X and Y coordinates
 
% Step 4: Create a 2D scatter plot with a heatmap color map
figure;
scatter(proj_coords(:, 2), proj_coords(:, 1), [], colors, 'filled');
xlim([-L/2, L/2])
 
% Apply color mapping using imagesc
colormap('jet');  % Choose a color map
shading interp;  % Interpolate colors between vertices
colorbar;  % Display color bar
 
xlabel('X');
ylabel('Time');
title('2D Projection with Heatmap');
 
% Add color bar
colorbar;

end
%%
function RHSvhat=rhs(dt,vhat)
% v should be a row vector
% RHSvhat = - e^{-tL}(1i*k*hat{(e^{tL}v)^2/2} 
N=size(vhat,2);
L = 32*pi;
k=-N/2 : (N/2 - 1);
freq =k.*(2*pi/L);
freq2 = freq.^2;
freq4 = freq2.^2;
e3=exp((freq2 - freq4)*dt); % in the Fourier space, uhat = e3.*vhat
em3=exp((-freq2 + freq4)*dt); % in the Fourier space, uhat = e3.*vhat
vhat1=vhat.*e3;          % e^{tL}v in the Fourier space 
v1=ifft(ifftshift(vhat1));      % exp(tL)v in the x-space
v2=0.5*v1.^2;          % [exp(tL)v]^2 in the x-space
RHSvhat=-em3.*(1i*freq).*fftshift(fft(v2)); % exp(-tL)[[(exp(tL)v)]_x] in the Fourier space
end
