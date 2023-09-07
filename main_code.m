% done by: Mohamed Magdy Abdelsatar 
%forward kinematics & Dynamics
% equation of motion of manipulators
% using symbolic toolbox and cell arrays:
data=fopen('datascara.txt');
scandata=textscan(data,'%s %s %s %s %s %s %s %s %s %s %s %s %f ','Delimiter',',');
fclose(data);
n =sym(scandata{1,13}(1,1)) % number of links
g = str2sym(scandata{1,12}(1,1)) %gravity matrix
double n;
M=cell(1,n);
iIi_all=cell(n,1);
irci_all=cell(n,1);
mi_all=cell(n,1);
qi_all=cell(n,1);
qi_dot_all=cell(n,1);
qi_double_dot_all=cell(n,1);
epseloni_all=cell(n,1);
for idx=1:n
    x='link';
    disp(x);
    disp(idx);    
    syms ai alphai thetai di mi qi qi_dot qi_double_dot
    syms iIi irci
    ai =str2sym(scandata{1}(idx));
    alphai =str2sym(scandata{2}(idx));
    di =str2sym(scandata{3}(idx));
    thetai =str2sym(scandata{4}(idx));
    syms a alpha d theta
    A=[ cos(theta)  -cos(alpha)*sin(theta)  sin(alpha)*sin(theta)   a*cos(theta);
        sin(theta)  cos(alpha)*cos(theta)   -sin(alpha)*cos(theta)  a*sin(theta);
        0           sin(alpha)              cos(alpha)              d;
        0           0                       0                       1;];
    % link i
    Ai=subs(A,{a,alpha,d,theta},{ai,alphai,di,thetai})
    M{1,idx}=Ai;
    %cell array
    mi=str2sym(scandata{5}(idx));
    mi_all{idx,1}=mi;
    iIi=str2sym(scandata{6}(idx));
    iIi_all{idx,1}=iIi;
    irci= str2sym(scandata{7}(idx));
    irci_all{idx,1}=irci;
    qi=str2sym(scandata{8}(idx));
    qi_all{idx,1}=qi;
    qi_dot=str2sym(scandata{9}(idx));
    qi_dot_all{idx,1}=qi_dot;
    qi_double_dot=str2sym(scandata{10}(idx));
    qi_double_dot_all{idx,1}=qi_double_dot;
    epseloni=sym(scandata{11}(idx));
    double epseloni
    epseloni_all{idx,1}=epseloni;
end
T = M{1} ;
A_all=cell(1,n);
A_all{1}=T;
R_0_i=cell(1,n);
R_0_i{1}=A_all{1}(1:3,1:3);
Zi_all=cell(1,n);
Zi_all{1}=A_all{1}(1:3,3);
for i = 2:numel(M)
    T = T * M{1,i}; %or: M{i}
    A_all{1,i}=T;
    R_0_i{1,i}=A_all{1,i}(1:3,1:3);
    Zi_all{1,i}=A_all{1,i}(1:3,3);
end
R_0_i
A_all
I_0_i_all=cell(1,n);
for x=1:n
    I_0_i=(R_0_i{1,x})*(iIi_all{x,1})*(transpose(R_0_i{1,x}));
    I_0_i_all{1,x}=I_0_i;
end
I_0_i_all
%link jacobian
J_v=cell(1,n);
J_w=cell(1,n);
z_alt=sym([0;0;1]);
z_zeros=sym([0;0;0]);
for y=1:n
    J_v_i=sym(cell(3,n));
    J_w_i=sym(cell(3,n));
    for w=1:n
        if w > y
            J_v_i(:,w)=z_zeros;
            J_w_i(:,w)=z_zeros;
        else
            if epseloni_all{w,1} == 0
                if w-1 == 0
                    J_v_i(:,w)=z_alt;
                    J_w_i(:,w)=z_zeros;
                else
                    J_v_i(:,w)=Zi_all{1,w-1};
                    J_w_i(:,w)=z_zeros;
                end
            else
                if w-1 == 0
                   J_w_i(:,w)=z_alt; 
                   J_v_i(:,w)=cross(z_alt,simplify(A_all{1,y}(1:3,:)*irci_all{y,1}));
                else
                   J_w_i(:,w)=Zi_all{1,w-1};
                   J_v_i(:,w)=cross(Zi_all{1,w-1},simplify((A_all{1,y}(1:3,:)*irci_all{y,1})-(A_all{1,w-1}(1:3,4))));
                end
            end
        end
    end
    J_v{1,y}=J_v_i;
    J_w{1,y}=J_w_i;
end
J_v
J_w
%Mass matrix
Mass=cell(n,n);
for s=1:n
    Mass_i=simplify((mi_all{s,1})*(transpose(J_v{1,s}))*(J_v{1,s}))+((transpose(J_w{1,s}))*(I_0_i_all{1,s})*(J_w{1,s}));
    if s==1
        Mass=Mass_i;
    else
        Mass=Mass+Mass_i;
    end
end
Mass
Mass_expanded=expand(Mass)
%velocity matrix
B=cell(1,n);
for f=1:n
    B_i=simplify((diff(Mass,qi_all{f,1}))*qi_dot_all);
    B{1,f}=B_i;
end
B
B_mat=[B{:}];
C=B-(0.5.*transpose(B_mat));
V=simplify(C*qi_dot_all);
V
%gravity matrix
G=cell(n,1);
for k=1:n
    for l=1:n
        G_i=simplify((-1).*(mi_all{l,1})*(transpose(g))*(J_v{1,l}(:,k)));
        if l==1
            G{k,1}=G_i;
        else
            G{k,1}=G{k,1}+G_i;
        end
    end
end
G
Mqdd=simplify(Mass*qi_double_dot_all)
%motion equation
tau=Mqdd+V+G
% tau=simplify(tau)