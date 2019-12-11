function [pop] = simulate_decay(init_pop, half_lives, t)
% Author: Ishbel Jamieson

% This function simulates the decay of D elements given the initial population,
% half lives, and time grid. 

% Input/output vectors in this function are sorted with the first element 
% decaying to the second element, then to the third element, and so on 
% until it reaches the last element - which is stable.

% This simulated decay is later compared against a simple analytic model
% for any observed deviations from theory. 

% Input:
% * init_pop: D-elements vector of the initial population of each
%             element, sorted such that the stable element is last.
% * half_lives: (D-1)-elements vector of half lives of each element.
% * t: N-elements vector of the time grid where the population is returned.

% Output:
% * pop: (N x D) matrix that contains the population of each element in 
%         every time step.

% Constraints:
% * The half life of an element can be very small, so giving the
%   probability to be (lambda*timestep) becomes an unsuitable model.
%   As such I have replaced this with (1 - exp(-lambda*timestep, which is
%   a more accurate model for elements with a large lambda value.
%
% Example use:
% >> init_pop = [1000, 0, 0, 0, 0];
% >> half_lives = [2.5e5, 8e4, 1.62e5, 4/365];
% >> t = linspace(0, 1e6, 1000);
% >> pop = simulate_decay(init_pop, half_lives, t);
%
% Example input:
% >> pop = simulate_decay([1000, 0, 0, 0, 0], [2.5e5, 8e4,1.62e5, 3.2e4], linspace(0,1e6,1000));


tstep = t(2) - t(1);
% Using a dummy variable zz, I will continue to redefine it throughout
% the function.
[zz,n] = size(t);
[zz,d] = size(half_lives);
% Assuming the last element is stable, it won't have a half life and
% therefore there will be d+1 elements to consider.
pop = zeros(n,1+d);
pop(1,:) = init_pop;
% Creating a probability column array and filling it by looping over d -
% therefore producing a more compact code. Note again that there is no given
% probability for the final stable element.
prob = zeros(d,1);
for p = 1:d
    %prob(p) = log(2)/half_lives(p)*tstep
    prob(p) = (1- exp(-(log(2)/half_lives(p)*tstep)));
end
% Again, the individual element loops are replace with one more general loop
% such that the code is more concise. A constraining line of code
% (line 53) was also required to ensure the total number of particles
% remained constant - without this the final element was under no boundary
% conditions and could therefore fluctuate between 1 and 0.
for time = 1:n-1
    pop(time+1,d+1) = pop(time,d+1);
    for elm = 1:d
        %Specifying the population of a particular element at a particular
        %time.
       for popu = 1:pop(time,elm)
           if prob(elm) >= rand
              pop(time+1,elm+1) = pop(time+1,elm+1) + 1;
           else
               pop(time+1,elm) = pop(time+1,elm) + 1;
           end 
       end
    end
end

% Creating an array for decay constant and filling it with a loop.
lambdas = zeros(d,1);
for i = 1:d
    lambdas(i) = log(2)/half_lives(i);
end

% Creating a matrix for the coefficients - ie the initial populations. Which
% we will then post-multiply by the expontential array, still in vector form.
% This will produce another matrix and due to the order of matrix
% multiplication - the i^th row of the new matrix is the expression for
% the population of the i^th element.

A = zeros(d-1);
for k = 1:d
    % The diagonal elements will be the term for the i^th element in the
    % expression for the the i^th element's population.
    A(k,k) = init_pop(k);
    for i = 1:k-1
        A(k,i) = A(k-1,i) * lambdas(k-1)/(lambdas(k)-lambdas(i));
        A(k,k) = A(k,k) - A(k,i);
    end
end

% Creating an analytical population matrix and setting the intial
% conditions.
anpop = zeros(n,d+1);
total_pop = sum(init_pop);
% Using an anonymous function (@(x)) inorder to take the exponent of each
% element in my array as I don't have MATLAB's Symbolab extension to deal
% with temporarily undefined functions/variables.
for time = 1:n
    expo = arrayfun(@(x) exp(x), -t(time)*lambdas);
    prod = A*expo;
    % Since 'prod' produces a d by n matrix, it must be transposed such that
    % any given row will now be the population of each element at a specific
    % time.
    anpop(time,1:d) = prod';
    % Separate line of code for the final stable element as it was not
    % considered in 'prod' (since I wanted to loop the array and the final
    % element doesn't follow the same pattern).
    anpop(time,d+1) = total_pop - sum(anpop(time,1:d));
end

% Comparing the simulation of a radioactive decay (pop) with the predicted
% trend line of any radioactive decay (anpop). (Whilst both adhere to our
% boundary conditions).
for r = 1:5
    plot(t,pop(:,r));
    hold on
    plot(t,anpop(:,r), '.');
    title('Comparison of Analytic and Simulated Radioactive Decay')
    xlabel('Time, years')
    ylabel('Number of Element Type in Population')
    % The element names are just illustrative of ones you could compare.
    legend('Uranium','Thorium','Radon','Radium','Lead')
    hold on
end
% From the plot, we can see the order in which the elements decay by
% comparing the time at which their maximum population occurs. We can also
% see how the simulation deviates from the theory on a microscopic level,
% though appears to be consistent on a macroscopic level.
end

