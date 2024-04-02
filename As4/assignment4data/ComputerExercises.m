%% CE 1
Clear all

A = imread("a.jpg");
B = imread("b.jpg");

figure()
subplot(1,2,1), imshow(A)
title('Image A')
subplot(1,2,2), imshow(B)
title('Image B')


[fA, dA] = vl_sift( single(rgb2gray(A))); 
[fB, dB] = vl_sift( single(rgb2gray(B)));
matches = vl_ubcmatch(dA,dB); % 204 matches
xA = fA(1:2,matches(1,:)); 
xB = fB(1:2,matches(2,:));

sample_size = 4;

bestH = [];
least_error = inf;
n_inliners = 0;
for steps = 1:10
    error = 0;
    Hi = DLT(xA, xB, sample_size);
    n = size(xA,2);
    p = 0; % Temporary inliner-counter
    for point=1:n
        sample = [xA(:,point);1];
        estimate = Hi*sample;
        estimate = estimate/estimate(end, end);
        point_error = sqrt((estimate(1)-xB(1,point))^2 + (estimate(2)-xB(2,point))^2);
        if point_error < 6
            p = p + 1;
        end
        error = error + point_error;
    end
    error = error/n;
    if p > n_inliners
        bestH = Hi;
        least_error = error; % least error 56.42 pixels
        n_inliners = p; %139 inliners found
    end
end

Htform = projective2d(bestH');

Rout = imref2d(size(A),[-200 800],[-400 600]);

[Atransf] = imwarp(A,Htform,"OutputView",Rout);

Idtform = projective2d(eye(3));
[Btransf] = imwarp(B,Idtform,"OutputView",Rout);

AB = Btransf;
AB(Btransf < Atransf) = Atransf(Btransf < Atransf);

figure()
imagesc(Rout.XWorldLimits ,Rout.YWorldLimits ,AB);

%% CE 2
close all
clc
load compEx2data.mat;
X1 = x{1};
X2 = x{2};

X1 = [X1; ones(1,size(X2,2))];
X2 = [X2; ones(1,size(X2,2))];

norm_X1 = inv(K)*X1;
norm_X2 = inv(K)*X2;

sample_size = 5;
best_inliner = -inf;
inliners = [];
best_E = [];
for iter=1:100
    ridx = randperm(size(X1, 2), sample_size);
    
    Es = fivepoint_solver(norm_X1(:,ridx), norm_X2(:,ridx));
    for e = 1:size(Es, 2)
        E = Es{1,e};
        F = (K^-1)'*E*K^-1;
        %De-normalize E to compute the distance in pixels
        %[average_error, p, good_points] = line_distance(X1, X2, pflat(F'*X1), pflat(F'*X2), 5);
        [average_error, inls, n_in] = line_distance2(pflat(X1), pflat(X2), pflat(F'*X2), pflat(F*X1), 5);
        if best_inliner < n_in
            best_E = E;
            best_inliner = n_in;
            inliners = inls;
            %best_points = good_points;
        end
    end
end

clearvars n_in average_error Es ridx F E e iter

% FROM ASSIGNMENT 3
[U,S,V] = svd(best_E);
if det(U*V') < 0
    V = -V;
end    

P1 = [diag([1,1,1]), zeros(3, 1)];
W = [0,-1,0;1,0,0;0,0,1];
u3 = U(:,3);

P21 = [U*W*transpose(V), u3];
P22 = [U*W*transpose(V), -u3];
P23 = [U*transpose(W)*transpose(V), u3];
P24 = [U*transpose(W)*transpose(V), -u3];

nPoints = size(norm_X1, 2);
% X3D = zeros(4, nPoints); % Initialize a matrix to hold 3D points

for i = 1:nPoints
    A1 = [norm_X1(1,i) * P1(3,:) - P1(1,:);
         norm_X1(2,i) * P1(3,:) - P1(2,:);
         norm_X2(1,i) * P21(3,:) - P21(1,:);
         norm_X2(2,i) * P21(3,:) - P21(2,:)];
    A2 = [norm_X1(1,i) * P1(3,:) - P1(1,:);
         norm_X1(2,i) * P1(3,:) - P1(2,:);
         norm_X2(1,i) * P22(3,:) - P22(1,:);
         norm_X2(2,i) * P22(3,:) - P22(2,:)];
    A3 = [norm_X1(1,i) * P1(3,:) - P1(1,:);
         norm_X1(2,i) * P1(3,:) - P1(2,:);
         norm_X2(1,i) * P23(3,:) - P23(1,:);
         norm_X2(2,i) * P23(3,:) - P23(2,:)];
    A4 = [norm_X1(1,i) * P1(3,:) - P1(1,:);
         norm_X1(2,i) * P1(3,:) - P1(2,:);
         norm_X2(1,i) * P24(3,:) - P24(1,:);
         norm_X2(2,i) * P24(3,:) - P24(2,:)];
    [~, ~, V1] = svd(A1);
    [~, ~, V2] = svd(A2);
    [~, ~, V3] = svd(A3);
    [~, ~, V4] = svd(A4);
    x1(:, i) = V1(:, end); 
    x2(:, i) = V2(:, end); 
    x3(:, i) = V3(:, end); 
    x4(:, i) = V4(:, end);
end

clearvars nPoints i

clearvars A1 A2 A3 A4 V1 V2 V3 V4 
x1 = pflat(x1);
x2 = pflat(x2);
x3 = pflat(x3);
x4 = pflat(x4);

dep1 = calculateDepth(P21, x1);
dep1 = [dep1, calculateDepth(P1,x1)];
dep2 = calculateDepth(P22, x2);
dep2 = [dep2, calculateDepth(P1,x2)];
dep3 = calculateDepth(P23, x3);
dep3 = [dep3, calculateDepth(P1,x3)];
dep4 = calculateDepth(P24, x4);
dep4 = [dep4, calculateDepth(P1,x4)];

% Camera 3 has 1472 points infront of it.
[positiveDepth, id] = max([sum(dep1 > 0),sum(dep2 > 0),sum(dep3 > 0),sum(dep4> 0)]);
clearvars dep1 dep2 dep3 dep4

if id == 1
    P_best = P21;
    x_best = x1;
end
if id == 2
    P_best = P22;
    x_best = x2;
end
if id == 3
    P_best = P23;
    x_best = x3;
end
if id == 4
    P_best = P24;
    x_best = x4;
end
clearvars P21 P22 P23 P24 x1 x2 x3 x4 

P1 = K*P1;
P_best = K*P_best;
x_hat = pflat(P1*x_best);
x_best_hat = pflat(P_best*x_best);

figure()
plot3(x_best(1,inliners == 1),x_best(2,inliners == 1),x_best(3,inliners == 1), "b.")
hold on
plotcams({P1, P_best})
axis equal

P = {P1, P_best};
[err, res] = ComputeReprojectionError(P, x_best, x);
RMS = sqrt(err/size(res,2));

figure()
hist(res,100)

%% CE 3
gamma = 1*10^(-10); % I don't understand how I'm supposed to work with gamma

Pnew = P;
Unew = x_best(:,inliners == 1);
x{1} = x{1}(:,inliners == 1);
x{2} = x{2}(:,inliners == 1);

[err, res] = ComputeReprojectionError(Pnew, Unew, x);

figure
plot(0, sum(res), 'b+')
hold on

for i = 1:10
    [r, J] = LinearizeReprojErr(Pnew, Unew, x);
    deltav = -gamma*J'*r;
    [Pnew, Unew] = update_solution(deltav, Pnew, Unew);
    [err , res] = ComputeReprojectionError(Pnew, Unew, x);
    plot(i, sum(res), 'b+');
end
RMS_ML = sqrt(err/size(res,2))

%% CE 4
lambda = 0.7; % I don't understand how I'm supposed to work with gamma

Pnew = P;
Unew = x_best(:,inliners == 1);

[err, res] = ComputeReprojectionError(Pnew, Unew, x);

figure
plot(0, sum(res), 'b+')
hold on

for i = 1:10
    [r, J] = LinearizeReprojErr(Pnew, Unew, x);
    C = J'*J+lambda*speye(size(J,2)); 
    c = J'*r;
    deltav = -C\c;
    [Pnew, Unew] = update_solution(deltav, Pnew, Unew);
    [err , res] = ComputeReprojectionError(Pnew, Unew, x);
    plot(i, sum(res), 'b+');
end
RMS_LM = sqrt(err/size(res,2))

function depth = calculateDepth(P, X)
    A = P(:, 1:3);
    a3 = P(3, end);
    X = X/X(end);
   
    signDetA = sign(det(A));
    
    normA3 = norm(A(3, :));
    
    product = [A(3, :), a3]*X;
    
    depth = (signDetA / normA3) * product;
end

function N = normise(x)
    % Calculate the mean of the points
    m = mean(x, 2);

    % Calculate the standard deviation of the points
    s = std(x, 0, 2);

    % Construct the normalization matrix
    N = [1/s(1), 0, -m(1)/s(1); 
         0, 1/s(2), -m(2)/s(2); 
         0, 0, 1];
end

function H = DLT(a, b, sample_size)
    i = randperm(size(a, 2), sample_size);
    a_sample = a(:,i);
    b_sample = b(:,i);

    M = zeros(2*sample_size,9);
    
    for sample=1:sample_size
        % [-x, -y, 1, 0,0,0, x'x, x'y,x']
        mx = [-a_sample(1,sample), -a_sample(2,sample), -1, 0, 0, 0, b_sample(1,sample)*a_sample(1,sample), b_sample(1,sample)*a_sample(2,sample), b_sample(1,sample)];
        my = [0, 0, 0, -a_sample(1,sample), -a_sample(2,sample), -1, b_sample(2,sample)*a_sample(1,sample), b_sample(2,sample)*a_sample(2,sample), b_sample(2,sample)];

        M(2*sample-1,:) = mx;
        M(2*sample,:) = my;
    end
    [~, ~, V] = svd(M);
    h = V(:, end);
    H = [h(1:3)'; h(4:6)'; h(7:9)'];
end

function [d, inliers, inlier_points] = line_distance(x, F, threshold)
    d = 0;
    n = size(X2,2);
    inliers = 0;
    inlier_points = [];
    lines1 = lines1 ./ sqrt(repmat(lines1(1, :).^2 + lines1(2, :).^2 ,[3 1]));
    lines2 = lines2 ./ sqrt(repmat(lines2(1, :).^2 + lines2(2, :).^2 ,[3 1]));

    for i = 1:n
        x2 = pflat(x{2}(:,i));
        x1 = pflat(x{1}(:,i));
        l1 = pflat(lines1(:, i));
        l2 = pflat(lines2(:, i));
        dist1 = abs(l1' * x2) / sqrt(l1(1)^2 + l1(2)^2); 
        dist2 = abs(l2' * x1) / sqrt(l2(1)^2 + l2(2)^2);

        if (dist1 < threshold) && (dist2 < threshold)
            inliers = inliers + 1; % Count each inlier pair once
            inlier_points = [inlier_points, i];
        end
        d = d + dist1 + dist2; 
    end
    d = d / (2 * n); % Correct average distance calculation
end

function [d, inls, inliers] = line_distance2(X1, X2, lines1, lines2, threshold)
    n = size(X2,2);
    x1 = X1;
    x2 = X2;
    inliers = 0;
   
    lines1 = lines1 ./ sqrt(repmat(lines1(1, :).^2 + lines1(2, :).^2 ,[3 1]));
    lines2 = lines2 ./ sqrt(repmat(lines2(1, :).^2 + lines2(2, :).^2 ,[3 1]));

    d1 = abs(sum(lines1.*x1)); % d1 = distance(l1, x1);
    d2 = abs(sum(lines2.*x2)); % d2 = distance(l2, x2);
    
    inls = (d1 < threshold) & (d2 < threshold);
    
    inliers = sum(inls(:));
    d = (d1+d2)/(2*n);
end


