

function [CL,on,onCVI]=DCVIClustering(varargin)
D=varargin{1}; alpha=2;

if length(varargin) == 2
     alpha=varargin{2};
end

[result,D1,LPS,CP,CL]=DCVI(D,alpha);
[on,onCVI,cl]=TestCVI(result);
CL(LPS)=cl;
for ii=1:size(D,1)
    if(CL(ii)==0)
        CL(ii)=CL(CP(ii));
    end
end
end

function[on,onCVI,cl]=TestCVI(varargin)
result=varargin{1};
kk=size(result,2)+1;
on=1;
onCVI=1;
for ii=2:ceil(sqrt(kk))
    [avgDCVI,bins]=testWeight(ii,result);
    temp=min(onCVI,avgDCVI);
    if onCVI~=temp
        onCVI=temp;
        on=ii;
        cl=bins;
    end
end
end
function[avgDCVI,bins,sortmst,weight2]=testWeight(varargin)
k=varargin{1};
mst=varargin{2};
avgDCVI=0;
DCVI=zeros(1,k);
edge=mst(3,:);
[Y,I]=sort(edge,2,'descend');
sortmst=mst(1:2,I);
sortmst(3,:)=Y;
cutnumber=k-1;
weight=zeros(1,k);
weight2=zeros(1,k);
cutmst=sortmst;
wmst=sortmst;

G=graph(cutmst(1,:),cutmst(2,:));
for ii=1:cutnumber
    a=cutmst(1,ii);
    b=cutmst(2,ii);
    G=rmedge(G,a,b);
end

bins=conncomp(G);
maxedge=sortmst(:,1:cutnumber);
cutpoint=zeros(3,cutnumber);
wmst(:,1:cutnumber)=[];
for ii=1:2
    for jj=1:cutnumber
        cutpoint(ii,jj)=bins(maxedge(ii,jj));
        
    end
end
cutpoint(3,:)=maxedge(3,:);

sep(1:k)=cutpoint(3,1);


for ii=1:k
    for jj=1:cutnumber
        if(cutpoint(1,jj)==ii||cutpoint(2,jj)==ii)
            sep(ii)=min(cutpoint(3,jj),sep(ii));
        end
    end
    
%     disp(sep(ii));
end
% disp('*******************************');

for ii=1:k
    temp2=0;
    temp=0;
    count=0;
    for jj=1:size(wmst,2)
        if bins(wmst(1,jj))==ii
            
            temp2=max(temp2,wmst(3,jj));
            temp=temp+wmst(3,jj);
            count=count+1;
            %       disp(temp2);
        end
    end
    if((count-1)==0)
        weight(ii)=0;
        
    else
        weight(ii)=temp/(count-1);
    end
    weight2(ii)=temp2;
    
end
temp=0;
for ii=1:k
    if (weight2(ii)==0)
        
        temp=min(weight2(weight2~=0));
        weight2(ii)=temp;
        
    end
%     disp(weight2(ii));
    
end

for ii=1:k
    DCVI(ii)=weight2(ii)/sep(ii);
end
avgDCVI=sum(DCVI)/k;

end

function [result,D1,LPS,CP,CL]=DCVI(varargin)
D=varargin{1}; alpha=varargin{2};
[CL,LPS,CP,CPV] = findCore(D,alpha);% Find the density core
% showCPV(D,CPV,LPS,CL);
D1 = D(LPS,1:end);         % Core point
[G,result] = Kruskal(D1);
end
function[CL,LPS,CP,CPV]=findCore(varargin)
D=varargin{1};alpha=varargin{2};

% Perform "NaN (Natural Neighbor) Searching"
[Sup,NN,RNN,NNN,nb,A]=NaNSearching(D);
[CPV,CP,Pr1,Pr2,r1,rf,r2,Nei1,CL]=findCenter2(D,NN,NNN,A,nb,alpha);

[LPS,FLP,T2]=findDensityPeak(CP,D,r1,rf,Pr1,Pr2,nb,Sup,CL);

end % Find the density core

function[Sup,NN,RNN,NNN,nb,A]=NaNSearching(varargin)
D=varargin{1};
r=1;
nb=zeros(size(D,1),1);
C=cell(size(D,1),1);
NN=cell(size(D,1),1);% Initialize the KNN neighbors of each point
RNN=cell(size(D,1),1);% Initialize the RKNN neighbors of each point
NNN=cell(size(D,1),1);% Is the intersection of NN and RNN, that is, for each point
A=pdist2(D,D);
Numb1=0;
Numb2=0;
for ii=1:size(D,1)
    [sa,index]=sort(A(:,ii));
    C{ii}=[sa,index];
end
while(r<size(D,1))
    for kk=1:size(D,1)
        x=kk;
        y=C{x}(r+1,2);
        nb(y)= nb(y)+1;
        NN{x}=[NN{x},y];
        RNN{y}=[RNN{y},x];
    end
    Numb1=sum(nb==0);
    if Numb2~=Numb1
        Numb2=Numb1;
    else
        break;
    end
    r=r+1;
end
for jj=1:size(D,1)
    NNN{jj}=intersect(NN{jj},RNN{jj});
end
Sup=r;
end
function[CPV,CP,Pr1,Pr2,r1,rf,r2,Nei1,CL]=findCenter2(varargin)
D=varargin{1};
NN=varargin{2};
NNN=varargin{3};
A=varargin{4};
nb=varargin{5};
alpha=varargin{6};


r1=[];
r2=[];
rf=[];
Pr1=[];
Pr2=[];
CPV=zeros(size(D,1),1);% Used for noise point inspection

CP=[];


Nei1=cell(size(D,1),1);
Nei2=cell(size(D,1),1);

for kk=1:size(D,1)
    CL(kk)=0;
end
for ii=1:size(D,1)
    if ~isempty(NN{ii})
        %         r2(ii)=mean(pdist2(D(ii,:),D(NN{ii},:)));% Take the mean distance between the point p and the surrounding natural neighbors as the radius to find the convergent point
        r1(ii)=1*max(pdist2(D(ii,:),D(NN{ii},:)));
        % Take the maximum value of the distance between the point p and the surrounding natural neighbors as the radius for seeking the density core
        r2(ii)=max(pdist2(D(ii,:),D(NN{ii},:)));
        rf=r1*0.95;
        Nei1{ii}=find(A(:,ii)<r1(ii));
        Nei2{ii}=find(A(:,ii)<rf(ii));
        Pr1(ii)=size(Nei1{ii},1);
        Pr2(ii)=size(Nei2{ii},1);
    else
        r1(ii)=0;
        r2(ii)=0;
        rf=r1*0.95;
    end
end
%
% Mark the noise points and prevent them from being mistakenly selected into RCP because the dynamic radius is too large. When there are no noise points in the data set, it is not necessary
B=mean(r2)+alpha*std(r2);
for ii=1:size(D,1)
    if r2(ii)>B
        CL(ii)=-1;
    end
    if r2(ii)==0
        CL(ii)=-1;
    end
    if nb(ii)<2
        CL(ii)=-1;
    end
end

for jj=1:size(D,1)% Clear the abnormalities in the points within the scanning radius of each point
    Nei1{ii}(find(CL(Nei1{ii})==-1))=[];
end

for ii=1:size(D,1)
    if ~isempty(Nei1{ii})
        if CL(ii)~=-1
            [~,y]=min(pdist2(D(Nei1{ii},:),mean(D(Nei1{ii},:),1)));
            if (CPV(Nei1{ii}(y))==ii) % Special treatment for points without neighbors and two points that converge to each other
                CPV(ii)=ii;
            else
                CPV(ii)=Nei1{ii}(y);
            end
        else
            CPV(ii)=ii;
        end
    else
        CPV(ii)=ii;
    end
end

for ii=1:size(D,1)% find the convergent point for point pi CP(i)
    if CL(ii)~=-1
        CP(ii)=ii;
        while(CP(ii)~=CPV(CP(ii)))
%             disp(CP(ii));
            CP(ii)=CPV(CP(ii));
        end
    else
        CP(ii)=ii;
    end
end

end
function[LPS,FLP,T2]=findDensityPeak(varargin)
CP=varargin{1};
D=varargin{2};
r1=varargin{3};
rf=varargin{4};
Pr1=varargin{5};
Pr2=varargin{6};
nb=varargin{7};
Sup=varargin{8};
CL=varargin{9};
LPS=[];
FLP=[];
T2=[];
for ii=1:size(D,1)
    if(CL(ii)~=-1)
        if(CP(ii)==ii)
            LPS=[LPS,ii];
            % %             T2=size(D,2)*log(rf(ii)/r1(ii))+log(Pr1(ii));
            %             if nb(ii)< Sup/2
            %                   FLP=[FLP,ii];
            %              end
        end
    end
    
end
end
function[]=showCPV(varargin)
D=varargin{1};
CPV=varargin{2};
LPS=varargin{3};
CL=varargin{4};
figure(1)
plot(D(CL~=-1,1),D(CL~=-1,2),'g.','markersize',5);
hold on
for ii=1:size(D,1)
    if (CPV(ii)~=ii&&CL(ii)~=-1)
        x=[D(ii,1),D(CPV(ii),1)];
        y=[D(ii,2),D(CPV(ii),2)];
        plot(x,y,'-','color',[0.4 0.8 0.9],'markersize',20)
    end
end

plot(D(LPS,1),D(LPS,2),'r.','markersize',10);
end% Draw the data set
function [T,result]=Kruskal(varargin)
% Construct minimum spanning tree
% Kruskal algorithm
D=varargin{1};
%% Make connection graph (adjacency matrix)
a = zeros(size(D,1));
for i = 1:size(D,1)
    for j = 1:size(D,1)
        if i<j
            a(i,j)=sqrt((D(i,1)-D(j,1))*(D(i,1)-D(j,1))+(D(i,2)-D(j,2))*(D(i,2)-D(j,2)));
        end
    end
end

G = graph(a,'upper');
[T,pred] = minspantree(G,'Method','sparse');
result = [ T.Edges.EndNodes(:,1), T.Edges.EndNodes(:,2),T.Edges.Weight]';
end % Construct minimum spanning tree
function[]=MST(varargin)
D=varargin{1};
result=varargin{2};
figure(2)
plot(D(:,1),D(:,2),'r.','markersize',10);
hold on


for i=1:size(result,2)
    x=[D(result(1,i),1),D(result(2,i),1)];
    y=[D(result(1,i),2),D(result(2,i),2)];
    plot(x,y,'b-');
    hold on
    
end


end % Draw minimum spanning tree