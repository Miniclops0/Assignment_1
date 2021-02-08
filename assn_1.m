% Connor Warden
% 101078296
% Assignment 1

%%% QUESTION LOCATION %%% 
%  Q_1_1
%  line 56
%  Q_1_2
%  line 59-60
%  Q_1_3
%  line 64-247 (includes most of the assignment)
%  Q_2_1
%  line 98-99, 102-106
%  Q_2_2
%  line  39, 237-244 
%  Q_2_3
%  line - shown in graph, varies 
%  Q_2_4
%  not completed
%  Q_3_1
%  line 43-48 (uncomment to turn on), 182-231
%  Q_3_2
%  line 52-53 (to turn on) 214-229 for implementation
%  Q_3_3
%  line 252-255
%  Q_3_4
%  line 262-268 - unfortunately does not work
%%%
clear all
close all

temp = 300; % 300 kelvin, given
b_constant = 1.38064852e-23; % Boltzmann's constant
m_o = 9.11e-31;
m_n = 0.26*m_o; % effective mass of electron, given

width = 200e-9; % width of region = 200 nano meters, given
height = 100e-9; % height of region = 100 nano meters, given
scattering = 0; % set = 1 for scattering
% Q_3_1
% Define Boxes, Uncomment for Question 3
Boxes = {}; % box list
% Boxes{1}.x = [0.8 1.2]*1e-7; % each cell element contains a structure
% Boxes{1}.y = [0.6 1.0]*1e-7; % starts at 60 nano meters and continues to the top
% 
% 
% Boxes{2}.x = [0.8 1.2]*1e-7; % each cell element contains a structure
% Boxes{2}.y = [0.0 0.4]*1e-7; % starts at 0 and continues to 40 nano meters, gap of 20 nano meters between the two boxes.

% Define if the boxes are specular or diffusive, set the preference to
% true and other to false
spec = true;
diff = false;

%Q_1_1
v_t = sqrt(2*b_constant*temp/m_n); % thermal velocity equation

%Q_1_2
tau_mn = 0.2e-12; 
lambda = v_t * tau_mn; % mean free path, for scattering

%Q_1_3

num_elec = 1000; % number of electrons simulated, only the first 10 in simulation

% Random particle location within x/y plane
x = rand(num_elec, 1)*width; % random x position 
y = rand(num_elec, 1)*height; % random y position 

% Makes sure nothing generates within the boxes, comment out for Q_1 and
% Q_2
for e = 1:num_elec
    if (0.8e-7 < x(e)) && (x(e) < 1.2e-7) % checks if within the box
        if x(e) > 1e-7 % if greater than the x/2 then move to the right
            x(e) = x(e) + 0.4e-7;
        else
            x(e) = x(e) - 0.4e-7; % if less than the x/2 then move to the left
        end
    end
end

xp = x; % x will be changed so we need to store the values
yp = y; % y will be changed so we need to store the values

step = 100; % number of iterations
dt = 0.01e-12; % time step
c = dt/v_t;

% Setting X and Y velocities

%Uncomment the following four lines for Q_1
r = rand(num_elec,1);
ang(:,1) = r(:,1)*2*pi;
vx = v_t*cos(ang); % maintains the same velocity but creates a different
vy = v_t*sin(ang); % launch angle for each iteration

% Uncomment the following 2 lines for Q_2 and Q_2
% vx = sqrt(b_constant*temp/m_n)*randn(num_elec, 1); % create x velocity
% vy = sqrt(b_constant*temp/m_n)*randn(num_elec, 1); % create y velocity

% Q_2_1
figure(1)
histogram(sqrt(vx.^2 + vy.^2), 100) % Generates histogram
title('Particle Velocity')
xlabel('Speed (m/s)')
ylabel("Number of Particles")



% Q_1 and Q_2 temp plot time steps
itr_a = zeros(1, 1);
for i = 1:num_elec-1 % creates an array from 0 - step*dt
    itr = i*dt;
    itr_a = cat(1, itr_a, itr);
end


% Check temps before, here just to check
v_avg = sqrt(vx.*vx + vy.*vy);
temp_out = ((v_avg.^2)*m_n)/(2*b_constant); 

% Main loop, handles every iteration

for i = 1:step 
        
    dx = vx*dt; % change in x

    dy = vy*dt; % change in y

    x = xp + dx; % updates the new x position
    y = yp + dy; % updates the new y position

    v_avg = sqrt(vx.*vx + vy.*vy);  % total avg velocity

    
    new_temp = ((v_avg.^2)*m_n)/(2*b_constant); % new temp value
    avg_temp = mean(new_temp, 'all'); % avg temp, displayed on plot
    figure(2)
    plot(itr_a, new_temp)
    
    
    % If the boxes are full, draw the box in the first iteration of
    % the for loop. Box is an input to signify this within the plot_traj fn
    if (i == 1) && ~isempty(Boxes)
        box = 1; % draw box
    else
        box = 0; % don't draw
    end

    colour = hsv(10); % colourmap 
    
    % Plots electron path
    % colour = array of colours
    % xp, yp = original pos    
    % x, y = new pos        
    % bix = if 1, plot boxes    
    plot_traj(colour, xp, x, yp, y, box)   
    
    % Define left/right top/bottom boundaries
    % y cannot be LESS than 0 or GREATER than the total height
    ylt = y < 0; 
    ygt = y > height;
    
    % x cannot be LESS than 0 or GREATER than the total width
    xlt = x < 0;
    xgt = x > width;
    
    % if the top/bottom are reached then reflect, i.e. invert the velocity 
    vy(ylt) = abs(vy(ylt));
    vy(ygt) = -abs(vy(ygt));
    
    % if sides are reached the move to the other side
    % leave the velocity in this instance, just change position
    x(xlt) = x(xlt) + width;
    x(xgt) = x(xgt) - width;
    
%     xp = x;
%     yp = y; 

    % Q_3_1/2 Handles box boundaries for both specular and diffusive

    if ~isempty(Boxes) % returns true if Boxes is not empty
        for c = 1:length(Boxes) % will go through every box case
            
            % Left/Right Box Boundary, makes sure to cover every case where 
            % the electron cannot pass
            inside_l = (x > Boxes{c}.x(1)) & (x < Boxes{c}.x(2)) & (xp < Boxes{c}.x(1)) & (y > Boxes{c}.y(1)) & (y < Boxes{c}.y(2)) ; 
            inside_r = (x > Boxes{c}.x(1)) & (x < Boxes{c}.x(2)) & (xp > Boxes{c}.x(2)) & (y > Boxes{c}.y(1)) & (y < Boxes{c}.y(2)) ;
            
            % Top/bottom boundaries differ based on if it the top or bottom
            % box, so if statement was used to define this. This can be
            % changed if there are more boxes but in the case of just 2
            % this works
            if (c == 1)  
                inside_b = (x > Boxes{c}.x(1)) & (x < Boxes{c}.x(2)) & (yp < Boxes{c}.y(1)) & (y > Boxes{c}.y(1));
                inside_t = (x > Boxes{c}.x(1)) & (x < Boxes{c}.x(2)) & (yp > Boxes{c}.y(2));
            else
                inside_b = (x > Boxes{c}.x(1)) & (x < Boxes{c}.x(2)) & (yp < Boxes{c}.y(1));
                inside_t = (x > Boxes{c}.x(1)) & (x < Boxes{c}.x(2)) & (yp > Boxes{c}.y(2)) & (y < Boxes{c}.y(2));  
            end
                
            x(inside_l) = Boxes{c}.x(1); % position of boundary, left
            x(inside_r) = Boxes{c}.x(2); % position of boundary, right
            
            y(inside_b) = Boxes{c}.y(1); % position of boundary, bottom
            y(inside_t) = Boxes{c}.y(2); % position of boundary, top
            
           % rand_l = randn(length(num_elec, 1)) < Boxes{c}.x(1);
           % rand_r = randn(length(num_elec, 1)) < Boxes{c}.x(2);
        
            
           % This decides whether to made boundaries specular or diffusive
           % based on a decision made near the beginning of the code
            if spec == true
                % Specular will reflect
                vx(inside_l) = -vx(inside_l); % invert velocity, reflects left
                vx(inside_r) = -vx(inside_r); % invert velocity, reflects right

                vy(inside_b) = -vy(inside_b); % invert velocity, reflects down for box 1, up for 2
                vy(inside_t) = -vy(inside_t); % invert velocity, reflects down for box 2, up for 1
            elseif diff == true
                % Diffusive will give a new random velocity
                vx(inside_l) = -(sqrt(2*b_constant*temp/m_n))*abs(randn(length(num_elec), 1)); % reflects left
                vx(inside_r) = (sqrt(2*b_constant*temp/m_n))*abs(randn(length(num_elec), 1)); % reflects right

                vy(inside_b) = -(sqrt(2*b_constant*temp/m_n))*abs(randn(length(num_elec), 1)); % invert velocity, reflects down for box 1, up for 2
                vy(inside_t) = (sqrt(2*b_constant*temp/m_n))*abs(randn(length(num_elec), 1)); % invert velocity, reflects down for box 2, up for 1
                
            end
        end
    end
    
    xp = x; % updates xp so that new path can be made on next itr
    yp = y; % updates yp so that new path can be made on next itr
    
%     % scattering
    if scattering == 1
        p = 1 - exp(-dt/tau_mn); % Probability of scattering
        r = p > rand(num_elec,1); % scattering value
        scat_l = length(vx(r));
    
        vx(r) = sqrt((b_constant*temp)/m_n)*randn(scat_l,1); % new velocity
        vy(r) = sqrt((b_constant*temp)/m_n)*randn(scat_l,1); % new velocity
    end
    
     
end

figure(2) % This just adds the average temperature once everything is done
title(['Semiconductor Temp ', num2str(avg_temp), ' K'])

% Plots the Electron Density Map
figure(4)
hist3([x, y], 'CDataMode','auto','FaceColor','interp')
title(' Electron Density')
% xlim([0, 200e-9])
% ylim([0, 100e-9])
% xlabel('x (m)')
% ylabel('y (m)')

% Plots the temperature map
figure(5)
for i = 1:x
    v_2 = mean(vx.*vx + vy.*vy);
    temp_map(i,i) = v_2*m_n/(2*b_constant);
    surf(temp_map)
    hold on
end













