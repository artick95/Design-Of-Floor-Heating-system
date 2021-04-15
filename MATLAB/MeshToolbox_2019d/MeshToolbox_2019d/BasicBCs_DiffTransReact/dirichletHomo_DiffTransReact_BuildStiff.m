function [D,T,R,b] = dirichletHomo_DiffTransReact_BuildStiff(Me, f)
%Assemble the matrix D and the vector b of the Transport-Diffusion-Reaction
%problem with homogeneous B.C.s
%Input:
%   Me     :a Mesh2D object
%   f      :MATLAB function of (x,y) which returns the values of the
%           external source. Default: constant value=4
%
%Output:
%   D      :diffusion matrix
%   T      :transport matrix
%   R      :reaction matrix
%   b      :constant terms vector

if nargin<2
    f=@(x,y)4*ones(size(x));
end 
%for clarity, call some properties of Me with shorter names
V=Me.Triangles.Vertices;
Areas=Me.Triangles.Areas;
CenterOfMass=Me.Triangles.CenterOfMass;
Nodes=Me.Nodes;
Dof=Me.Nodes.Dof;
%number of internal nodes: we know that the N unknown nodes are numbered from
%1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
%(degrees of freedom)
numDof = max(Dof);

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
b = zeros(numDof,1);
row = zeros(Me.MatrixContributions,1);
col = zeros(Me.MatrixContributions,1);
r=zeros(Me.MatrixContributions,1);
t=zeros(Me.MatrixContributions,1);
d=zeros(Me.MatrixContributions,1);
pos=1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry
%evaluate the value of the coefficient in front of the Laplace operator
mu=Me.mu;
%evaluate the value of the coefficient in front of the divergence operator
%NB. beta is a two-columns matrix, indicating the speed in the x
%and y directions
beta=Me.beta;
rho=Me.rho;
%evaluate the value of the advection/reaction coefficient
sigma=Me.sigma;
%we evaluate the external force in the center of mass of each triangle
force = f(CenterOfMass.X,CenterOfMass.Y);
%main loop on each triangle        
for e=1:size(V,1)   
    Dx(1) = Nodes.X(V(e,3)) - Nodes.X(V(e,2));
    Dx(2) = Nodes.X(V(e,1)) - Nodes.X(V(e,3));
    Dx(3) = Nodes.X(V(e,2)) - Nodes.X(V(e,1));
    Dy(1) = Nodes.Y(V(e,3)) - Nodes.Y(V(e,2));
    Dy(2) = Nodes.Y(V(e,1)) - Nodes.Y(V(e,3));
    Dy(3) = Nodes.Y(V(e,2)) - Nodes.Y(V(e,1));
    
    %for each vertex of this triangle 
    for ni=1:3
        %look at the "unknown" numbering: if the node is positive, it
        %corresponds to a degree of freedom of the problem
        ii = Dof(V(e,ni));
        %is it unknown?
        if ii > 0 
            %yes it is! second loop on the vertices
            for nj=1:3
                jj = Dof(V(e,nj));
                %%is it unknown as well?
                if jj > 0
                    %Non sparse solution: 
                    %D(ii,jj)=D(ii,jj) + c*(Dy(i)*Dy(j)+Dx(i)*Dx(j))/(4.0*Area) ;
                    %T(ii,jj)=T(ii,jj) + (beta(1)*Dy(j)-beta(2)*Dx(j))*1/6;
                    %R(ii,jj)=R(ii,jj) + a*Area*(1+(i==j))/12;
                    row(pos)=ii;
                    col(pos)=jj;
                    t(pos)=(beta(e,1)*Dy(nj)-beta(e,2)*Dx(nj))*1/6*rho(e);
                    d(pos)=mu(e)*(Dy(ni)*Dy(nj)+Dx(ni)*Dx(nj))/(4.0*Areas(e));
                    r(pos)=sigma(e)*Areas(e)*(1+(ni==nj))/12;
                    pos=pos+1;                    
                end
            end
            %build the constant terms vector adding the external
            %contribution
            b(ii) = b(ii)+ Areas(e)*force(e)/3.0;
        end
    end
end
%assemble the stiffness matrix D from the
D=sparse(row, col, d, numDof, numDof);
T=sparse(row, col, t, numDof, numDof);
R=sparse(row, col, r, numDof, numDof);
