% Newton Divided Difference Interpolation

diary cos1

l=[1,-1,2,-2,3,-3,4,-4,14,-14,1000,-1000];
n = size(l,2);
for i=1:1:n
    fprintf('x = %d\n', l(i))
    fprintf('cos(x): %.7g\n', cos(l(i)))
    val = cos1([0,pi/4,pi/2], [1, 1/sqrt(2),0],l(i));
    fprintf('cos1(x): %.7e\n\n', val)
    fprintf('Error: %.7e\n\n', abs(val-cos(l(i))))
end

diary off

function val = cos1(x,y, x_o)
    n = size(x,1);
    n2 = size(x,2);
    if n < n2
        x = x';
        y = y';
        n = n2;
    end
    coefs = newtonDDsetup(x, y);
    val = coefs(1);
    val = val + coefs(2)*(x_o-x(1)) + coefs(3)*(x_o-x(1))*(x_o-x(2));

end


function newtonDD1(x, y)
    n = size(x,1);
    n2 = size(x,2);
    if n < n2
        x = x';
        y = y';
        n = n2;
    end
    
    coefs = newtonDDsetup(x, y);
    disp('x='); xt=x';disp(xt)
    disp('y='); yt=y';disp(yt)
    coefs    
    
    m = 10*n;
    minx = min(x);
    maxx = max(x);
    t = minx:(maxx-minx)/m:maxx;
    
    val = newtonEval(t, coefs, x);
    figure
    plot(t,val);
    hold on;
    for i=1:1:n
        plot(x(i),y(i),'k*');
    end
    xlabel('x'); ylabel('y');
    title('Newton DD interpolation');
    legend('Newton DD interpolant','Data points','Location','best');
end
    
    
% Evaluate divided difference interpolant
function value = newtonEval(t, coefs, x)
    n = size(x,1);
    value = coefs(n);
    for i=n-1:-1:1
        value = value .* (t - x(i)) + coefs(i);
    end
end


% Set up divided difference coefficients
function coefs = newtonDDsetup(x, y)
    n = size(x,1);

    % DD level 0
    for i=1:1:n
        coefs(i) = y(i);
    end

    % DD higher levels (bottom to top, overwrite lower entries as they are finished)
    for level=2:1:n
        for i=n:-1:level
            dx = x(i) - x(i-level+1);
            coefs(i) = ( coefs(i)-coefs(i-1) ) / dx;
        end
    end
end




