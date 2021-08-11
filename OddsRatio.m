classdef OddsRatio 
  properties
    a,b,c,d;
  end

  methods
    function obj=OddsRatio(a,b,c,d)
    % New odds ratio:  a/b vs c/d
      obj.a=a;
      obj.b=b;
      obj.c=c;
      obj.d=d;
    end
    
    function ratio=getRatio(obj)
      ratio=(obj.a/obj.b)/(obj.c/obj.d);
    end
    
    function [alpha,pr]=getAlphaBoot(obj,nboot)
      if nargin<2
        nboot=1000;
      end
      univ=false(obj.a+obj.b+obj.c+obj.d,1);
      univ(1:obj.a+obj.c)=true;
      n1=obj.a+obj.b;
      n2=obj.c+obj.d;
      pr=[];
      for i=1:nboot
        pick=univ(randperm(length(univ)));
        pr(i)=sum(pick(1:n1))/sum(~pick(1:n1))/(sum(pick(n1+1:end))/sum(~pick(n1+1:end)));
      end
      r=obj.getRatio();
      req=(obj.a+obj.c)/(obj.b+obj.d);
      if r>req
        alpha=mean(pr>=r)*2;
      else
        alpha=mean(pr<=r)*2;
      end
    end
    

    function alpha=getAlphaNorm(obj)
    % Calculate alpha using a normal approximation
      n1=obj.a+obj.b; p1=(obj.a+0.5)./(n1+1);v1=p1.*(1-p1)./n1;
      n2=obj.c+obj.d; p2=(obj.c+0.5)./(n2+1);v2=p2.*(1-p2)./n2;
      % Distribution of difference in p
      dp=abs(p1-p2);
      v=v1+v2;
      cum=normcdf(dp,0,sqrt(v));
      alpha=(1-cum)*2;
    end
    
    function [low,high]=getCI(obj,alpha,varargin)
    % Get confidence interval for (a/b) / (c/d)
      defaults=struct('method','fleiss');
      args=processargs(defaults,varargin);

      a=obj.a; b=obj.b; c=obj.c; d=obj.d;
      
      if strcmp(args.method,'fleiss')
        za2=norminv(alpha/2);
        m=a+c; n=b+d;
        s=a+b; f=c+d;
        N=m+n;

        for k=1:2
          if k==1
            v=obj.getRatio()/1.1;
            c=-1/2;
          else
            v=obj.getRatio()*1.1;
            c=1/2;
          end
          for i=1:5000
            if mod(i,100)==1
              fprintf('k=%d, i=%d, v=%f\n', k, i, v);
            end
            X=v*(m+s)+(n-s);
            Y=sqrt(X^2-4*m*s*v*(v-1));
            A=(X-Y)/(2*(v-1));
            B=s-A;
            C=m-A;
            D=f-m+A;
            W=1/A+1/B+1/C+1/D;
            F=(a-A+c)^2*W-za2^2;
            T=1/(2*(v-1)^2)*(Y-n-(v-1)/Y*(X*(m+s)-2*m*s*(2*v-1)));
            U=1/B^2+1/C^2+1/A^2+1/D^2;
            V=T*((a-A+c)^2*U-2*W*(a-A+c));
            v=v-F/V;
            assert(isfinite(v));
          end
          if k==1
            low=v;
          else
            high=v;
          end
        end
      else
        error('Unknown method: %s',args.method);
      end
    end

  end % methods
end % class