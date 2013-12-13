function savefigure(directory, filename)
    % SAVEFIGURE
    %   saveFigure(directory,filename) saves the figure as a png, jpg, eps,
    %   pdf and fig all at once ***
    %
    % ***: check current file: some saving formats may be commented out
    %
    
    directory = fixpath(directory);
    
    % if the directory doesn't exist, ask for it to be created
    dircontents = dir(directory);
    if L(dircontents)==0
        fprintf(['Directory: ' directory ' does not exist \n']);
        fprintf('Do you want to:\n');
        fprintf('   - 1: create it\n');
        fprintf('   - 2: reenter it\n');
        fprintf('   - 3: cancel?\n');
        aa = input('                  >> ');
        if aa==1
            mkdir(directory);
        elseif aa==2
            newdir = input('enter new directory:  ','s');
            savefigure(newdir,filename);
            return;
        elseif aa==3
            return;
        else
            fprintf('invalid input. try again\n\n\n');
            savefigure(directory,filename);
            return;
        end
    end
    
    print('-dpng','-r300',[directory filename '.png']);
    %print('-djpeg','-r125',[directory filename ' - PREVIEW.jpg']);
    %print('-depsc',[directory filename '.eps']);
    %print('-dpdf',[directory filename '.pdf']);
    hgsave([directory filename '.fig']);


end