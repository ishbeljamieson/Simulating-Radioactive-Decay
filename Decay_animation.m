% Author: Ishbel Jamieson

% This script visually simulates radioactive decay with an animation. Each
% of the three colours correspond to a different element, with the final
% element Z (represented by blue pixels) being the stable one.


% Example Input:
% Pixel depth = 20
% Decay rate of X = 1
% Decay rate of Y = 1
% Timestep = 0.25
% Endyear = 200


% Creating a 3D matrix that will form the image, the size of third
% dimension of 'M' indicates that it's an RGB image.
p = str2double(inputdlg({'Pixel depth:'}));
M = zeros(p,p,3);
input = inputdlg({'Decay rate of X:','Decay rate of Y:','Timestep:','End year:'}, 'ALL values must be in years');

% I will associate the following RGB colour codes to each element:
% 110 (Yellow) represents element X
% 010 (Green) represents element Y
% 001 (Blue) represents element Z

% Setting all initial pixels to be yellow, our original element X.
% element.
M(:,:,1)=1;
M(:,:,2)=1;
M(:,:,3)=0;
image(M)

%Assigning the inputed data to variables.
decay_x = str2double(input(1,:));
decay_y = str2double(input(2,:));
timestep = str2double(input(3,:));
endyear = str2double(input(4,:));



for n = 1:round(endyear/timestep)
% Looping over pixel depth to fill the rows and columns.
for r = 1:p
    for c = 1:p
        if M(r,c,1)==1 && M(r,c,2) == 1 && M(r,c,3) ==0
           x = rand;
           prob_decay_x = 1 - exp(-1*decay_x*timestep);
           if prob_decay_x > x
              M(r,c,1)=0;
           end
           
        elseif M(r,c,1) == 0 && M(r,c,2) == 1 && M(r,c,3) == 0
               x=rand;
               prob_decay_y = 1 - exp(-1*decay_y*timestep);
               if prob_decay_y > x
                  M(r,c,2) = 0;
                  M(r,c,3) = 1;
               end
        end
    end
end
pause(1)
image(M)
a = 0;
b = 0;
d = 0;
for r = 1:p
    for c = 1:p
        if M(r,c,1)==1 && M(r,c,2) == 1 && M(r,c,3) == 0
           a = a + 1;
        elseif M(r,c,1) == 0 && M(r,c,2) == 1 && M(r,c,3) == 0
               b = b + 1;
        else
            d = d + 1;
        end
    end
end

end
