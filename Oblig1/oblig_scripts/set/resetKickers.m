% SET ALL KICKERS TO 0 (inactive)
function [] = resetKickers()

    nKickers = 8;
    kicks0 = zeros(1,nKickers);
    setKickers(kicks0, kicks0)
    
end

