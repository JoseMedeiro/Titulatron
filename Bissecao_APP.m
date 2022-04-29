       
function D = Bissecao_APP(DATA_OLD, DATA_NEW)
    if DATA_NEW(size(DATA_NEW,1),1) < 0
        if DATA_OLD(size(DATA_OLD,1),1) < 0
            D = [DATA_NEW , DATA_OLD(:,2)];
        else
            D = [DATA_OLD(:,1) , DATA_NEW];
        end
    else
        if DATA_OLD(size(DATA_OLD,1),1) > 0
            D = [DATA_NEW , DATA_OLD(:,2)];
        else
            D = [DATA_OLD(:,1) , DATA_NEW];
        end
    end
end