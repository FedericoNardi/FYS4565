% SET KICKERS
function [] = setKickers( xKicks, yKicks )

    % define kicker file used by MADX
    kicker_file = 'kickers.tfs';

    % save kicks
    fileID = fopen(kicker_file,'w');
    fprintf(fileID, '@ NAME %%07s "KICKERS"\n');
    fprintf(fileID, '@ TYPE %%04s "USER"\n');
    fprintf(fileID, '* XKICK YKICK\n');
    fprintf(fileID, '$ %%le %%le\n');
    for i = 1:numel(xKicks)
        fprintf(fileID, ' %4.2e %4.2e \n', [xKicks(i), yKicks(i)]);
    end
    fclose(fileID);

end
