classdef DisturbanceLinearSystem < LinearSystem

    properties (SetAccess = private)
        W % convex set of distrubance
        Z % disturbance invariant set
    end

    methods (Access = public)

        function obj = DisturbanceLinearSystem(A, B, Q, R, W)
            obj = obj@LinearSystem(A, B, Q, R);
            obj.W = W;
            Nmax=25;
            obj.Z = obj.compute_mrpi_set(Nmax);
        end

        function x_new = propagate(obj, x, u)
            w = obj.pick_random_disturbance();
            x_new = propagate@LinearSystem(obj, x, u) + w;
        end

    end

    methods (Access = public)

        function w = pick_random_disturbance(obj)
            % pick disturbance form uniform distribution
            verts = obj.W.V;
            b_max = max(verts)';
            b_min = min(verts)';
            w = rand(obj.nx, 1) .* (b_max - b_min) + b_min; 
            % generate random until it will be inside of W
%             while true
%                 w = rand(obj.nx, 1) .* (b_max - b_min) + b_min; 
%                 if obj.W.contains(w)
%                     break
%                 end
%             end
        end

        function Fs = compute_mrpi_set(obj, Nmax)
%             % Computes an invariant approximation of the minimal robust positively
%             % invariant set for 
%             % x^{+} = Ax + w with w \in W
%             % according to Algorithm 1 in 'Invariant approximations of
%             % the minimal robust positively invariant set' by Rakovic et al. 
%             % Requires a matrix A, a Polytope W, and a tolerance 'epsilon'. 
            AkW = cell(Nmax + 1, 1);
            alphaW = cell(Nmax + 1, 1);
            alphas = NaN(Nmax + 1, 1);
            Ak=eye(2);
            verts = (obj.W.V)';
            for n = 1:(Nmax + 1)
                AkW{n} = Ak*verts;
                normW = norm(verts,inf);
                normKW=norm(obj.K*verts,inf);
                normAkW = norm(AkW{n},inf);
                normKAkW = norm(obj.K*AkW{n},inf);
                
                alphas(n) = max(normAkW/normW, normKAkW/normKW);
                alphaW{n} = alphas(n)*verts;
                Ak=obj.Ak*Ak;
                PAkw{n}=Polyhedron((AkW{n})');  
                if n>5
                    if max(max(abs(AkW{n}-AkW{n-1})))<0.001
                        break
                    end    
                end
            end
            
            Fs=PAkw{1};
            for i=2:n
                Fs=Fs+PAkw{i};
            end 
            Fs=Fs*(1/(1-alphas(n)));
         end

    end
end
