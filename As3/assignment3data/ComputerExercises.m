%% Computer Exercise 1
clear
close all
load compEx1data.mat

pic1 = imread("Kronan1.JPG");
pic2 = imread("Kronan2.JPG");

% xhatFx = 0
% Normalization matrix
N1 = normize(x{1}(1:2, :));
N2 = normize(x{2}(1:2, :));

%N1 = diag([1,1,1])
%N2 = N1

% Normailzed points
xn1 = N1*x{1};
xn2 = N2*x{2};

M = zeros(size(xn1,2), 9);

for i=1:size(xn1,2)
    M(i, :) = reshape(xn2(:, i) * xn1(:, i)', 1, []);
end

[~, ~, V] = svd(M);
v = V(:,end); %both are small
n = norm(M*v, 2); %both are small

% Finding the normalized F
%Fn = reshape(v, 3, 3)';
Fn = reshape(v, [3 3]);
[Uf, Sf, Tf] = svd(Fn);
Sf(3, 3) = 0; 
Fn = Uf * Sf * Tf';
plot(diag(transpose(xn2)*Fn*xn1));

% Going back to unnormalized
F = N2' * Fn * N1;
F = F./F(end,end);

det(F); % 2.8387e-18
constraint = transpose(x{2}(:,1))*F*x{1}(:,1); % Close to zero

% EPIPOLAR LINES
l = F*x{1};
l = l./sqrt(repmat(l(1,:).^2 + l(2,:).^2,[3 1]));
figure()
dist = histogram(abs(sum(l.*x{2})),100);

meandist = mean(abs(sum(l.*x{2}))); % with normalization = 23.92, without = 34.52
num_points = 20;
i = randi([1, size(xn1, 2)], 1, num_points); % 20 random points
points = x{2}(:,i);
lines = l(:,i);

figure();
imagesc(pic2);
hold on
rital(lines);
plot(points(1,:), points(2,:), 'y+', MarkerSize=30)

%% Computer Exercise 2
e2 = null(F');
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
P2 = [e2x * F, e2]; % Camera 2 in normalized form
P1 = [diag([1,1,1]), zeros(3, 1)];

P2n = N2*P2;
P1n = N1*P1;

nPoints = size(xn1, 2);
X = zeros(4, nPoints); % Initialize a matrix to hold 3D points

for i = 1:nPoints
   A = [xn1(1,i) * P1n(3,:) - P1n(1,:);
        xn1(2,i) * P1n(3,:) - P1n(2,:);
        xn2(1,i) * P2n(3,:) - P2n(1,:);
        xn2(2,i) * P2n(3,:) - P2n(2,:)];
   [~, ~, V] = svd(A);
   X(:, i) = V(1:4, end);
end


x1_h = pflat(P1n*X);
x2_h = pflat(P2n*X);
x1_h = inv(N1)*x1_h;
x2_h = inv(N2)*x2_h;

% P1 = inv(N1)*P1n;
% P2 = inv(N2)*P2n;
% 
% X = pflat(X);
% 
% x1_hat = pflat(P1*X);
% x2_hat = pflat(P2*X);

figure()
imagesc(pic1);
hold on
plot(x{1}(1,:), x{1}(2,:), "b.", MarkerSize=10);
plot(x1_h(1,:), x1_h(2,:), "r+", MarkerSize=8);

figure()
imagesc(pic2);
hold on
plot(x{2}(1,:), x{2}(2,:), "b.", MarkerSize=10);
plot(x2_h(1,:), x2_h(2,:), "r+", MarkerSize=8);

figure()
X = pflat(X);
plot3(X(1,:), X(2,:), X(3,:), '.')
title("3D-Plot")
hold on

%% Computer Exercise 3
load compEx3data.mat
xn1 = K^(-1)*x{1};
xn2 = K^(-1)*x{2};

% Same as 
M = zeros(size(xn1,2), 9);
for i=1:size(xn1,2)
    place = xn2(:,i)*xn1(:,i)';
    M(i,:) = place(:)';
end

[~, ~, V] = svd(M);
v = V(:,end); %both are small
n = norm(M*v); %both are small
Eapprox = reshape(v, 3, 3);
[U,~,V] = svd(Eapprox); % Two singular values and one zero
if det(U*transpose(V))>0
    E = U*diag([1 1 0])*transpose(V); 
else
    V = -V;
    E = U*diag([1 1 0])*transpose(V); 
end
constraint = transpose(xn1(:,1))*E*xn2(:,1); % Close to zero
E = E./E(3,3);

% e2 = null(E');
% e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];

%F = e2x*E;
%F = K*F./F(end,end);
F = transpose(K^(-1))*E*K^(-1)./F(end,end);
l = F*x{1};
l = l ./ sqrt(repmat(l(1, :).^2 + l(2, :).^2 ,[3 1]));
i = randi([1, size(xn1, 2)], 1, 20); % 20 random points
points = x{2}(:,i);
lines = l(:,i);

figure();
imagesc(pic2);
hold on
rital(lines);
plot(points(1,:), points(2,:), 'r+', MarkerSize=30)

figure()
dist = histogram(abs(sum(l.*x{2})),100);
meandist = mean(abs(sum(l.*x{2}))); % 2.803

%% Computer Exercise 4
[U,S,V] = svd(E);
P1 = [diag([1,1,1]), zeros(3, 1)];
W = [0,-1,0;1,0,0;0,0,1];
%U = [(1/sqrt(2)),(-1/sqrt(2)),0;(1/sqrt(2)),(1/sqrt(2)),0;0,0,1];
u3 = U(:,3);

P21 = [U*W*transpose(V), u3];
P22 = [U*W*transpose(V), -u3];
P23 = [U*transpose(W)*transpose(V), u3];
P24 = [U*transpose(W)*transpose(V), -u3];

nPoints = size(xn1, 2);
X = zeros(4, nPoints); % Initialize a matrix to hold 3D points

for i = 1:nPoints
    A1 = [xn1(1,i) * P1(3,:) - P1(1,:);
         xn1(2,i) * P1(3,:) - P1(2,:);
         xn2(1,i) * P21(3,:) - P21(1,:);
         xn2(2,i) * P21(3,:) - P21(2,:)];
    A2 = [xn1(1,i) * P1(3,:) - P1(1,:);
         xn1(2,i) * P1(3,:) - P1(2,:);
         xn2(1,i) * P22(3,:) - P22(1,:);
         xn2(2,i) * P22(3,:) - P22(2,:)];
    A3 = [xn1(1,i) * P1(3,:) - P1(1,:);
         xn1(2,i) * P1(3,:) - P1(2,:);
         xn2(1,i) * P23(3,:) - P23(1,:);
         xn2(2,i) * P23(3,:) - P23(2,:)];
    A4 = [xn1(1,i) * P1(3,:) - P1(1,:);
         xn1(2,i) * P1(3,:) - P1(2,:);
         xn2(1,i) * P24(3,:) - P24(1,:);
         xn2(2,i) * P24(3,:) - P24(2,:)];
    [~, ~, V1] = svd(A1);
    [~, ~, V2] = svd(A2);
    [~, ~, V3] = svd(A3);
    [~, ~, V4] = svd(A4);
    X1(:, i) = V1(:, end); 
    X2(:, i) = V2(:, end); 
    X3(:, i) = V3(:, end); 
    X4(:, i) = V4(:, end);
end
X1 = pflat(X1);
X2 = pflat(X2);
X3 = pflat(X3);
X4 = pflat(X4);

dep11 = calculateDepth(P21, X1);
dep12 = calculateDepth(P1, X1);
dep21 = calculateDepth(P22, X2);
dep22 = calculateDepth(P1, X2);
dep31 = calculateDepth(P23, X3);
dep32 = calculateDepth(P1, X3);
dep41 = calculateDepth(P24, X4);
dep42 = calculateDepth(P1, X4);

% Camera 2 has 2008 points infront of the cameras.
[positiveDepth, id] = max([sum(dep11 > 0 & dep12 >0),sum(dep21 > 0 & dep22 >0),sum(dep31 > 0 & dep32 >0),sum(dep41 > 0 & dep42 >0)]);

P1 = K*P1;
P21 = K*P21;
P22 = K*P22;
P23 = K*P23;
P24 = K*P24;
x_hat = pflat(P1*X1);
x1_hat = pflat(P21*X1);
x2_hat = pflat(P22*X2);
x3_hat = pflat(P23*X3);
x4_hat = pflat(P24*X4);

figure()
imagesc(pic1);
hold on
plot(x{1}(1,:), x{1}(2,:), "b.");
plot(x_hat(1,:), x_hat(2,:), "r+");

figure()
imagesc(pic2);
hold on
plot(x{2}(1,:), x{2}(2,:), "b.");
plot(x1_hat(1,:), x1_hat(2,:), "r+");

figure()
imagesc(pic2);
hold on
plot(x{2}(1,:), x{2}(2,:), "b.");
plot(x2_hat(1,:), x2_hat(2,:), "r+");

figure()
imagesc(pic2);
hold on
plot(x{2}(1,:), x{2}(2,:), "b.");
plot(x3_hat(1,:), x3_hat(2,:), "r+");

figure()
imagesc(pic2);
hold on
plot(x{2}(1,:), x{2}(2,:), "b.");
plot(x4_hat(1,:), x4_hat(2,:), "r+");


%% 3D plot
figure()
plot3(X2(1,:),X2(2,:),X2(3,:), "b.")
hold on
plotcams({P1, P22})
axis equal

function depth = calculateDepth(P, X)
    A = P(:, 1:3);
    a3 = P(3, end);
    X = X/X(end);
   
    signDetA = sign(det(A));
    
    normA3 = norm(A(3, :));
    
    product = [A(3, :), a3]*X;
    
    depth = (signDetA / normA3) * product;
end

function N = normize(x)
    m = mean(x, 2);

    s = std(x, 0, 2);

    N = [(1/s(1))   0       -(1/s(1))*m(1);
         0      (1/s(2))    -(1/s(2))*m(2);
         0      0       1];
end