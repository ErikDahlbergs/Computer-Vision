%% Exercise 1
% part 1
t_x = [0,0,2;0,0,0;-2,0,0];
A = [1,1,0; 0,2,0; 0,0,1];

F = t_x*A;

% part 2
x = [1;1;1];
l = F*x

%% Exercise 2
F = [0,0,2;0,0,-2;-2,2,-2];
% Part 1
C_1 = [0;0;0;1];
P_1 = [1,0,0,0;0,1,0,0; 0,0,1,0];

C_2 = [-1;-1;0;1];
P_2 = [1,1,1,2;0,2,0,2; 0,0,1,0];

e1 = P_1*C_2
e2 = P_2*C_1;

% Part 2

true1 = transpose(e2) * F
true2 = F*e1 

% Part 3
dets = det(F) % är 0

%% Exercise 4
F = [0,1,1;1,0,0;0,1,1];
%[U,S,V] = svd(F);
%e2 = V(:,end)./V(end,end);
e2 = [-1, 0, 1];
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
P2 = [e2x*F, e2'];
P2 = [-1 0 0 -1; 0 2 2 0; -1 0 0 1];
P1 = [diag([1,1,1]), zeros(3, 1)];
X = cell(3,1);
X{1} = [1; 2; 3];
X{2} = [3; 2; 1];
X{3} = [1; 0; 1];
x_h = cell(3,2);
constraint = zeros(3,1);
for i = 1:3
    x1 = P1 * [X{i}; 1];
    x2 = P2 * [X{i}; 1];
    x_h{i,1} = x1;
    x_h{i,2} = x2;
end

for i = 1:3
    % The points are already in homogeneous coordinates
    x1 = x_h{i,1};
    x2 = x_h{i, 2};

    % Calculate the epipolar constraint
    constraint(i) = x2' * F * x1;

    % Display the result
    fprintf('Constraint for point pair %d is: %f\n', i, constraint(i));
end


%% Exercise 6
% Part 1
U = [(1/sqrt(2)),(-1/sqrt(2)),0;(1/sqrt(2)),(1/sqrt(2)),0;0,0,1];
V_T = transpose([1,0,0;0,0,(-1);0,1,0]);

dets = det(U*V_T); % = 1 

% Part 2
diagonal = diag([1,1,0]);
E = U*diagonal*V_T;
x1 = [0;0;1];
x2 = [1;1;1];

verify = transpose(x2)*E*x1

% Part 3
W = [0,-1,0;1,0,0;0,0,1];

B1 = U*W*V_T;
B2 = U*transpose(W)*V_T;
u3 = [0 ; 0 ; 1];

S1 = [B1, u3];
S2 = [B1, -u3];
S3 = [B2, u3];
S4 = [B2, -u3];

% Part 4
% A3_1 = S1(3,:);
% A3_2 = S2(3,:);
% A3_3 = S3(3,:);
% A3_4 = S4(3,:);

X12= [-1/(sqrt(2));-1/(sqrt(2));-1/(sqrt(2));1];
X34 = [1/(sqrt(2));1/(sqrt(2));1/(sqrt(2));1];

X1= [0;0;1;-1/(sqrt(2))];
X2= [0;0;1;1/(sqrt(2))];
X3= [0;0;1;1/(sqrt(2))];
X4= [0;0;1;-1/(sqrt(2))];

depth1 = calculateDepth(S1, X1);
depth2 = calculateDepth(S2, X2);
depth3 = calculateDepth(S3, X3);
depth4 = calculateDepth(S4, X4);

P1 = [eye(3,3) zeros(3,1)];

d1_1 = calculateDepth(P1, X1)
d1_2 = calculateDepth(P1, X2)
d1_3 = calculateDepth(P1, X3)
d1_4 = calculateDepth(P1, X4)

function depth = calculateDepth(P, X)
    A = P(:, 1:3);
    a3 = P(3, end);
    X = X/X(end);
   
    signDetA = sign(det(A));
    
    normA3 = norm(A(3, :));
    
    product = [A(3, :), a3]*X;
    
    depth = (signDetA / normA3) * product;
end





