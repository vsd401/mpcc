

function x0 = unWrapX0(x0)

    if ( x0(3) > pi)
        x0(3) = x0(3) - 2*pi;
    end
    if( x0(3) <= -pi)
       x0(3) = x0(3) + 2*pi;
    end
end