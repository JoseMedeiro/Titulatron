
function D = DIST_REL(A, B)
    if (A+B~=0)
        D = abs(DIST_ABS(A,B)/(A + B)*2);
    else
        if A
            D = abs(DIST_ABS(A,B)/A);
        else
            D = 0;
        end
    end
end