function [op, nodes] = mk_nbrs_of_digraph_del_add(G0,A,pa_limit,hprior,numEdge,numParent)
% Generate all possible single edge modification(deletion/addition) for a given network
% INPUTS:
%	G0:	adj matrix s.t. G0(i,j)=1 iff i->j in graph
%	A: 	path cont matrix
%   pa_limit:   the number of maximum parents
%   hprior:    A matrix indicating allowed edges (1=allowed, 0=banned)
%   numEdge:    the number of edges in current graph
%   numParent:    the parents number of each nodes in a current graph
% OUTPUTS:
%	op:		operation for each edge (0=delete,2=add)
%   nodes: source node and target node of the edges can be modified
%
% The Original code from bnt library
% https://code.google.com/p/bnt/
% Code covered by the GNU GPL v2 License
%
% The original script has been modified to only perform a single edge 
% deletion and addition, and increase speed, and incoporate hprior.
% Modification made on Nov. 2012 by:
%
% Shinya Tasaki, Ph.D.


%% SINGLE EDGE DELETIONS

% all edges are deletable
[I,J] = find(G0); % I(k), J(k) is the k'th edge

%% SINGLE EDGE ADDITIONS

% SML: previously Kevin had that any addition was legal if A(i,j)=0
% however, you can not add i->j  if j is a descendent of i.
% Thus, we create all possible additions in Gbar and then
% subtract the descendants of each edge as possible parents
% This means the potential parents of i (i.e. Gbar(:,i))
% can not also be descendants if i i.e. (A(:,i)) which is accomplished
% by subtracting (Gbar-A == 1 iff Gbar=1 & A=0)


GbarL=(G0-A')==0;
GbarL(:,numParent== pa_limit)=0;
[IbarL, JbarL] = find(GbarL&hprior);  % I(k), J(k) is the k'th legal edge to add
EbarL = length(IbarL);

%% put together
nodes = [I J;
   IbarL JbarL];
op = [zeros(numEdge,1);ones(EbarL,1)+1];
