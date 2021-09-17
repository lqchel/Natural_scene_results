clear all
try
    
    %% 1. Collect subject ID number and initials
    %Screen('Preference', 'SkipSyncTests', 1);
    
    cd('C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos')
    %% input subject ID and number
    subnumber=input('Enter Subject Number (01-15): \n','s');
    subj.initials = input('Enter subject initials:\n','s'); % 'AL'
    SesNumber = input('Enter Session Number (01/02): \n','s');
    BlockNumber = input('Enter Block Number (01/02): \n','s');
    %% load image and patches locations
    scene_dir = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\shinji images'; %% images for scrambled masks
    Incong_dir = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\squareimage\incongruent cropped'; %% incongruent image folder
    Cong_dir= 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\squareimage\congruent cropped'; %% congruent image folder
    testPatch_dir = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\patches'; %%  folder with both congruent and incongruent patches
    shinjiPatch_dir = 'C:\Users\Chelsea\Documents\honours\research project\Experiment\New photos\nishimoto patch'; %% null patches
    
    % Add auxiliary files folder to path
    addpath('C:\Users\Chelsea\OneDrive\Documents\honours\research project\Experiment\New photos\aux_files');
    
    %% create filename
    save_path = ['Expresult_' subj.initials '_' subnumber '_' SesNumber '_' BlockNumber '.mat'];
    
    
    %% Load pilot trial list and specify variables
    
    ExpList = importdata(['Triallist_', num2str(subnumber),'_', num2str(SesNumber),num2str(BlockNumber) '.mat']);
    PracticeList = importdata('practicelist.mat');
    if SesNumber == "01" && BlockNumber == "01"
        TrialList = [PracticeList; ExpList];
    else
        TrialList = ExpList;
    end
    %transform the trial number index in Trial List into numbers
    %for later trial list indexing
    for tnum = 1:length(TrialList)
        TrialNumberList(tnum,:) = str2num(TrialList(tnum,2));
    end
    %% INITIAL DEFINITIONS:
    %number of practice trials
    num_prac = 2;
    
    % number of experiment only trials
    num_exp_trials= 20;
    
    %number of total trials
    if SesNumber == "01" && BlockNumber == "01"
        num_trials= num_exp_trials+num_prac;
    else 
        num_trials = num_exp_trials;
    end
    
    block1 = 5;
    block2 = 10;
    block3 = 15;
    
    %% 3. Initiate Psychtoolbox Screen
    screen_size = 'full';
    % The practice_parameters function will now be called to "open" our window
    % using Psychtoolbox. Don't worry too much about this one, it prepares our
    % experiment for use with this screen and allows us to call and assign
    % things to it using 'Exp'
    
    Exp = parameters(screen_size);
    %% SOME FINAL DEFINITIONS
    
    % create duration of presentation
    % DEFINE LENGTH OF EACH SOA as NUMBER OF FRAMES.
    par.fr=Screen('FrameRate',Exp.Cfg.win);
    %e.g. frame= 1/refreshrate
    
    % for mac:
    % par.duration4=0.133333334;
    par.duration4 = 16/par.fr; %~133ms
    %CHOOSE DURATION FOR EXPERIMENT
    par.currentDuration=par.duration4;
    
    par.num_trials=num_trials;
    par.num_prac_trials=num_prac;
    par.num_expOnly_trials=num_exp_trials;
    num_masks=5;
    par.num_masks=num_masks;
    % Fixation cross presentation time (seconds)
    fixation_time = 0.5; % 500ms... might need to change this.
    par.fixation_time=fixation_time;
    % Time each mask is presented (seconds)
    mask_time = 0.06; %60ms
    par.mask_time=mask_time;
    
    %wait time
    wait_time_bw_screens=0.25;
    par.wait_time_bw_screens=wait_time_bw_screens;
    wait_time_responseWheel=0.08;
    par.wait_time_responseWheel=wait_time_responseWheel;
    
    %present or absent
    present=-1;
    absent=1;
    yes=1;
    no=-1;
    
    %% FIXATION CROSS SETTINGS
    
    %Set colour, width, length etc.
    CrossColour = 0;  %255 = white
    CrossL = 15;
    CrossW = 3;
    
    %Set start and end points of lines
    Lines = [-CrossL, 0; CrossL, 0; 0, -CrossL; 0, CrossL];
    CrossLines = Lines';
    
    
    % Anything starting with 'Exp' was defined in the 'parameters'
    % function.
    % We define what we want to happen to our screen then 'flip' the result so
    % these three lines tell us that we want to fill the screen with the
    % colour black, present some text in the middle of the screen, and then
    % 'flip' the result to the experiment window
    %% presentation of title
    
    intImage = imread('instruction.jpg');
    Instruction_Tex = Screen('MakeTexture',Exp.Cfg.win, intImage);
    Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.black);
    Screen('DrawTexture', Exp.Cfg.win, Instruction_Tex);
    Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
    
    % ...And this little while-loop waits for a mouse-click or button-press before we
    % continue
    while (1)
        [~,~,buttons] = GetMouse(Exp.Cfg.win);
        if buttons(1) || KbCheck
            break;
        end
    end
    
    WaitSecs(wait_time_bw_screens);
    
    
    %% 4. Show first trial
    
    for trial = 1:num_trials
        %% load test image
        exp_trial = trial - num_prac; % remove first 4 practice trials from all experiment images
        
        %find the current trial list in the master TrialList,
        %by comparing the trial number list with current trial number (i.e. trial)
        %which returns a logical array and index the target trial list location in the master list
        %and we use the logical array to extract the current trial list from the master list
        Im_ID = trial;

        TrialIndex = TrialNumberList == Im_ID;
        CurrentTrial = TrialList(TrialIndex,:);
        
        %based on the image feature (congruent/incongruent) specified in the trial list, load the corresponding image
        %load images for practice and exp trials
        
        if contains(CurrentTrial(1,3), "In")
            if Im_ID < 10
                temp = imread([Incong_dir, '\SquareIncongruent_00', num2str(Im_ID), '.jpg']);
            elseif Im_ID >= 10 && Im_ID<100
                temp = imread([Incong_dir,'\SquareIncongruent_0', num2str(Im_ID),'.jpg']);
            else
                temp = imread([Incong_dir,'\SquareIncongruent_', num2str(Im_ID),'.jpg']);
            end
            
        else
            if Im_ID < 10
                temp = imread([Cong_dir, '\SquareCongruent_00', num2str(Im_ID), '.jpg']);
            elseif Im_ID >= 10 && Im_ID<100
                temp = imread([Cong_dir,'\SquareCongruent_0', num2str(Im_ID),'.jpg']);
            else
                temp = imread([Cong_dir,'\SquareCongruent_', num2str(Im_ID),'.jpg']);
            end
        end
        temp=imresize(temp,[880 880]); %adjust image size
        
        
        
        
        
        %% image mask details setup
        %if trial > num_prac
        for maskCounter = 1:num_masks %change to length of variable
            mask_index = randi([1 500], 1, 1);
            filePattern1 = fullfile(scene_dir,'*.jpg');
            theFiles1 = dir(filePattern1);
            filename1 = theFiles1(mask_index).name;
            fullname1 = fullfile(scene_dir,filename1);
            rndmask = imread(fullname1);
            [mask,NewOrder(maskCounter,:)] = scramble(rndmask, rndmask); %take random image and make mask
            Mask_Tex(maskCounter) = Screen('MakeTexture',Exp.Cfg.win,mask);
        end
        
        %end
        
        %% define duration of presentation for image:
        
        flip_time = par.currentDuration;
        
        
        if trial > 1 && ~(trial == num_prac+1 || exp_trial == block1 || exp_trial == block2 || exp_trial == block3)
            Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.black);
            DrawFormattedText(Exp.Cfg.win,'Click mouse button when \n\n ready for the next trial.','center','center',Exp.Cfg.Color.white);
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            
            while (1)
                [~,~,buttons] = GetMouse(Exp.Cfg.win);
                if buttons(1) || KbCheck
                    break;
                end
            end
            WaitSecs(wait_time_bw_screens);
            
        elseif trial == num_prac+1
            Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.black);
            DrawFormattedText(Exp.Cfg.win,'PRACTICE END \n\n Press any key on keyboard to Start Actual Experiment','center','center',Exp.Cfg.Color.white);
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            while (1)
                [keyIsDown, secs, keycode] = KbCheck();
                kp = find(keycode);
                if ~isempty(kp)
                    break;
                end
            end
            WaitSecs(wait_time_bw_screens);
            
            
        elseif exp_trial == block1 || exp_trial == block2 || exp_trial == block3
            Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.black);
            DrawFormattedText(Exp.Cfg.win,['You have completed ' num2str(exp_trial) ' trials. \n\n Click mouse button when you are ready to continue.'],'center','center',Exp.Cfg.Color.white);
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            
            while (1)
                [~,~,buttons] = GetMouse(Exp.Cfg.win);
                if buttons(1) || KbCheck
                    break;
                end
            end
            WaitSecs(wait_time_bw_screens);
        end
        
        
        
        % First, let's hide the cursor
        HideCursor;
        
        % The code below will 'make' textures for the image and trial textures based upon those defined by the
        % 'create_trials' script
        
        Image_Tex = Screen('MakeTexture',Exp.Cfg.win, temp);
        
        
        %% FIXATION CROSS
        
        % Now, let's colour the screen gray
        Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
        Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
        
        % ... And show the fixation cross for 'fixation_time' amount of seconds
        time_remaining = fixation_time;
        targetSecs = GetSecs;
        
        % This while-loop will 'flip' until time_remaining is equal to zero
        while time_remaining > 0
            
            Screen('DrawLines', Exp.Cfg.win, CrossLines, CrossW, CrossColour,...
                [Exp.Cfg.xCentre, Exp.Cfg.yCentre]);
            
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            
            time_elapsed = GetSecs - targetSecs;
            time_remaining = fixation_time - time_elapsed;
            
        end
        
        %% changed by NT 17 June 12
        Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
        %%
        
        % Prepare while-loop timing
        time_remaining = flip_time; % Defined above
        startSecs = GetSecs; % Checks the current computer clock
        
        %% FLIP IMAGE TO SCREEN
        % changed by NT 17 June 12
        %    while time_remaining > 0
        
        isDrawn = 0 ;
        currentMask = 1;
        
        while GetSecs - startSecs < flip_time - 0.0042 % uncomment for linux: 1/par.fr / 2
            %        Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
            if isDrawn == 0
                Screen('DrawTexture', Exp.Cfg.win, Image_Tex);
                Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
                isDrawn = 1;
                
                % prepare the first mask here
                %scaScreen('DrawTexture', Exp.Cfg.win, Mask_Tex(currentMask));
                Screen('DrawTexture', Exp.Cfg.win, Mask_Tex(currentMask),[0 0 880 880]);
            end
            
            %        time_elapsed = GetSecs - startSecs; % Check time
            %        time_remaining = flip_time - time_elapsed; % Remove from time_remaining
            
        end
        
        %% MASKS x 5
        currentMask = 1;
        isTimedMaskStart = 0;
        % while currentMask <= num_masks
        for currentMask = 1:num_masks
            
            isMaskDrawn = 0;
            EachMaskStartTime = GetSecs ;
            %        while time_remaining_masks > 0
            while GetSecs - EachMaskStartTime < mask_time
                % time_remaining_masks > 0
                % Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
                
                if isMaskDrawn == 0
                    if currentMask ~= 1
                        Screen('DrawTexture', Exp.Cfg.win, Mask_Tex(currentMask), [0 0 880 880]);
                    end
                    Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
                    
                    if currentMask == 1 &&  isTimedMaskStart == 0
                        startSecs_masks_actual = GetSecs; % Checks the current computer clock
                        isTimedMaskStart = 1;
                    end
                    isMaskDrawn = 1;
                end
                
                % time_elapsed = GetSecs - startSecs_masks; % Check time
                % time_remaining_masks = mask_time - time_elapsed; % Remove from time_remaining
                
                % time_remaining_masks
                
            end
            
            %        currentMask = currentMask +1
        end
        %% Blank screen for 80ms
        Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
        Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
        
        WaitSecs(wait_time_responseWheel);
        clear rndmask
        
        %% PATCH PRESENTATION
        
        for pnum = 1:length(CurrentTrial) % pnum = patch number
            %% image info to be saved
            
            if trial == 1
                TestedPatchNum = pnum;
                
            else
                TestedPatchNum = TestedPatchNum + 1;
            end
            CurrentTrialNumber(TestedPatchNum,:) = trial;
            
            if CurrentTrial(1,3)== ['Congruent']
                ImageFeature(TestedPatchNum,:) = 0;
            else
                ImageFeature(TestedPatchNum, :) = 1;
            end
            
            %% fixation cross before patch presentation
            Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            
            % ... And show the fixation cross for 'fixation_time' amount of seconds
            time_remaining = fixation_time;
            targetSecs = GetSecs;
            
            % This while-loop will 'flip' until time_remaining is equal to zero
            while time_remaining > 0
                
                Screen('DrawLines', Exp.Cfg.win, CrossLines, CrossW, CrossColour,...
                    [Exp.Cfg.xCentre, Exp.Cfg.yCentre]);
                
                Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
                
                time_elapsed = GetSecs - targetSecs;
                time_remaining = fixation_time - time_elapsed;
                
            end
            
            
            
            %% load patch
            PatchID = char(CurrentTrial(pnum,4));
            
            if contains(PatchID, '.jpg')
                loadpatch = imread([testPatch_dir, '\',PatchID]);
                patch = imresize(loadpatch, 4);
            else
                filepattern = fullfile(shinjiPatch_dir,'*.jpg');
                theFiles = dir(filepattern);
                filename = theFiles(str2num(PatchID)).name;
                patchname = fullfile(shinjiPatch_dir, filename);
                loadpatch = imread(patchname);
                patch = imresize(loadpatch, 4);
            end
            %% patch mask setup
        for maskCounter = 1:num_masks %change to length of variable
            mask_index = randi([1 500], 1, 1);
            filePattern1 = fullfile(scene_dir,'*.jpg');
            theFiles1 = dir(filePattern1);
            filename1 = theFiles1(mask_index).name;
            fullname1 = fullfile(scene_dir,filename1);
            rndmask = imread(fullname1);
            [mask,NewOrder(maskCounter,:)] = scramble(rndmask, rndmask); %take random image and make mask
            PatchMask_Tex(maskCounter) = Screen('MakeTexture',Exp.Cfg.win,mask);
        end
            
            %% store patch source and features for saving responses
            if contains(PatchID, 'in')
                PatchSource(TestedPatchNum,1) = 3; %patch from incongruent image, coded as 3
            elseif contains(PatchID, 'cong')
                PatchSource(TestedPatchNum,1) = 2; %patch from congruent image, coded as 2
            else
                PatchSource(TestedPatchNum,1) = 1; %patch from nishimoto's image, coded as 1
            end
            
            if contains(PatchID, '_a')
                PatchFeature(TestedPatchNum,1) = 2; %critical object absent
            elseif contains(PatchID, '_p')
                PatchFeature(TestedPatchNum,1) = 3; %critical object present
            else
                PatchFeature(TestedPatchNum,1) = 1; %patch from nishimoto's image
            end
            
            QuestionNumber(TestedPatchNum, 1)= pnum;
            %% generate and set patch presentation locations
            %set 9 possible locations for patch presentation
            location_5 = CenterRect([0 0 880/3 880/3], Exp.Cfg.windowRect); %set centre location first, then draw it to other locations by distance on x- and y- axis
            location_1 = OffsetRect(location_5,-880/3, -880/3);
            location_2 = OffsetRect(location_5, 0, -880/3);
            location_3 = OffsetRect(location_5, 880/3, -880/3);
            location_4 = OffsetRect(location_5,-880/3, 0);
            location_6 = OffsetRect(location_5, 880/3, 0);
            location_7 = OffsetRect(location_5,-880/3, 880/3);
            location_8 = OffsetRect(location_5, 0, 880/3);
            location_9 = OffsetRect(location_5, 880/3, 880/3);
            
            %define presentation location for current patch
            if contains(PatchID,string([num2str(Im_ID),'_1_'])) %for patches from congruent/incongruent images, we want to present them in the same location as they are in the picture
                CurrentLocation = location_1;
                Patch_location(TestedPatchNum,1) = 1;
            elseif contains(PatchID,string([num2str(Im_ID),'_2_']))
                CurrentLocation = location_2;
                Patch_location(TestedPatchNum,1) = 2;
            elseif contains(PatchID,string([num2str(Im_ID),'_3_']))
                CurrentLocation = location_3;
                Patch_location(TestedPatchNum,1) = 3;
            elseif contains(PatchID,string([num2str(Im_ID),'_4_']))
                CurrentLocation = location_4;
                Patch_location(TestedPatchNum,1) = 4;
            elseif contains(PatchID,string([num2str(Im_ID),'_5_']))
                CurrentLocation = location_5;
                Patch_location(TestedPatchNum,1) = 5;
            elseif contains(PatchID,string([num2str(Im_ID),'_6_']))
                CurrentLocation = location_6;
                Patch_location(TestedPatchNum,1) = 6;
            elseif contains(PatchID,string([num2str(Im_ID),'_7_']))
                CurrentLocation = location_7;
                Patch_location(TestedPatchNum,1) = 7;
            elseif contains(PatchID,string([num2str(Im_ID),'_8_']))
                CurrentLocation = location_8;
                Patch_location(TestedPatchNum,1) = 8;
            elseif contains(PatchID,string([num2str(Im_ID),'_9_']))
                CurrentLocation = location_9;
                Patch_location(TestedPatchNum,1) = 9;
            else
                location_set = {location_1 location_2 location_3 location_4 location_5 location_6 location_7 location_8 location_9}; %for Nishimoto patch, we randomly select a presentation location from the location set
                random_locnum = randi([1,9], 1, 1);
                CurrentLocation = cell2mat(location_set(random_locnum));
                Patch_location(TestedPatchNum,1) = random_locnum;
            end
            
            
            
            %% present patch on corresponding location; location
            Patch_Tex = Screen('MakeTexture',Exp.Cfg.win, patch);
            
            isPresented = 0 ;
            PatchSecs = GetSecs;
            currentMask = 1;
            while GetSecs - PatchSecs < flip_time - 0.0042 % uncomment for linux: 1/par.fr / 2
                %        Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
                %present patch in corresponding location(location as showed in patch filename)
                if isPresented == 0
                    Screen('DrawTexture', Exp.Cfg.win, Patch_Tex, [0 0 880/3 880/3], CurrentLocation);
                    Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
                    isPresented = 1;
                    
                    %present mask
                    Screen('DrawTexture', Exp.Cfg.win, PatchMask_Tex(currentMask), [0 0 880/3 880/3], CurrentLocation);
                    % prepare the first mask here
                    %scaScreen('DrawTexture', Exp.Cfg.win, Mask_Tex(currentMask));
                    
                end
                
                %        time_elapsed = GetSecs - startSecs; % Check time
                %        time_remaining = flip_time - time_elapsed; % Remove from time_remaining
                
            end
            
            %% MASKS x 5
            currentMask = 1;
            isTimedMaskStart = 0;
            % while currentMask <= num_masks
            for currentMask = 1:num_masks
                %        if currentMask == 1
                %            startSecs_masks = GetSecs; % Checks the current computer clock
                %        end
                %        time_remaining_masks = mask_time; % Defined above
                isMaskDrawn = 0;
                EachMaskStartTime = GetSecs ;
                %        while time_remaining_masks > 0
                while GetSecs - EachMaskStartTime < mask_time
                    % time_remaining_masks > 0
                    % Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
                    
                    if isMaskDrawn == 0
                        if currentMask ~= 1
                            Screen('DrawTexture', Exp.Cfg.win, PatchMask_Tex(currentMask),[0 0 880/3 880/3], CurrentLocation);
                        end
                        Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
                        
                        if currentMask == 1 &&  isTimedMaskStart == 0
                            startSecs_masks_actual = GetSecs; % Checks the current computer clock
                            isTimedMaskStart = 1;
                        end
                        isMaskDrawn = 1;
                    end
                    
                    % time_elapsed = GetSecs - startSecs_masks; % Check time
                    % time_remaining_masks = mask_time - time_elapsed; % Remove from time_remaining
                    
                    % time_remaining_masks
                    
                end
                
                %        currentMask = currentMask +1
            end
            
            %% Blank screen for 80ms
            Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            
            WaitSecs(wait_time_responseWheel);
            
            %end
            
            
            
            %% DRAW RESPONSE SCREEN
            % At this point, all of our image textures have been shown so let's
            % present the response wheel so we can register a decision
            
            % Firstly, let's show the cursor again
            ShowCursor;
            
            % Now, we 'flip' the response wheel to our screen and wait a moment to
            % ensure people don't accidentally click (this stuff happens very
            % quickly)
            if trial > num_prac
                [Exp] = response_screen(Exp,'Do you think the patch presented to you belongs to the image?','  ', 'Yes','No');
            else
                [Exp] = response_screen(Exp,'Do you think the patch presented to you belongs to the image?','  ','Yes','No');
            end
            timeStartResponse=GetSecs;
            WaitSecs(wait_time_bw_screens);
            
            %% REGISTER A RESPONSE
            % We set up a while-loop that will break once we register a click
            clicks = 0;
            
            while clicks == 0
                
                [x,y] = getMouseResponse(); % This waits until a click is made
                
                if checkQuit %if participant press q to quit experiment, break while loop
                   break
                end
                % This stuff is a little confusing but it checks whether a click
                % went inside one of the boxes that are built into the confidence
                % wheel
                for m = 1 : size(Exp.polyL, 1)
                    idxs_left(m) = inpolygon(x,y,squeeze(Exp.polyL(m,1,:)),squeeze(Exp.polyL(m,2,:)));
                    
                    idxs_right(m) = inpolygon(x,y,squeeze(Exp.polyR(m,1,:)),squeeze(Exp.polyR(m,2,:)));
                end
                
                idx_pos_left = find(idxs_left == 1);
                idx_pos_right = find(idxs_right == 1);
                
                % Left boxes click == ANSWER = NO
                if length(idx_pos_left) == 1  %~isempty(idx_pos_left)
                    keyid = no; %answered no
                    keyid2 = idx_pos_left;
                    
                    clicks = 1;
                    
                    % Paint selected box red
                    Screen('FillPoly', Exp.Cfg.win, [255 0 0], squeeze(Exp.polyL(idx_pos_left,:,:))',1);
                    for wait = 1:10
                        Screen('Flip', Exp.Cfg.win,  [], Exp.Cfg.AuxBuffers);
                    end
                    
                end
                %%ANSWER = YES
                if length(idx_pos_right) == 1 %~isempty(idx_pos_right)
                    keyid = yes; %answered no
                    keyid2 = idx_pos_right;
                    
                    clicks= 1;
                    
                    % Paint selected box blue
                    Screen('FillPoly', Exp.Cfg.win, [0 0 255], squeeze(Exp.polyR(idx_pos_right,:,:))',1);
                    for wait = 1:10
                        Screen('Flip', Exp.Cfg.win,  [], Exp.Cfg.AuxBuffers);
                    end
                    
                end
                
            end
            
            timeEndResponse=GetSecs;
            
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck(); %if participants press q, break patch testing loop
            if keyCode(KbName('q'))
                
                break
                
            end
 
            %% 6. Save responses to .mat file
            % 1 is present, 2 is absent
            % 1 is left/no, -1 is right/yes
            %             TR(exp_trial).keyid(currentScreen) = keyid; %save yes or no
            
            % SIGNAL = 'presence'
            % DECISION or RESPONSE = keyid
            % Absolute CONFIDENCE = keyid2
            Response(TestedPatchNum,:)= keyid;
            Confidence(TestedPatchNum,:) = keyid2;
            %response = keyid .* keyid2
            
            % Measure Response Time from when the response screen is presented until subject makes a click.
            ResponseTime(TestedPatchNum, :) =timeEndResponse-timeStartResponse;
            
            
            Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            
            WaitSecs(wait_time_responseWheel);
            
           
            
            
        end   %end patch testing loop
        
     %% Check whether q is pressed; if pressed, break trial loop to abort experiment & save results  
        
         [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
            if keyCode(KbName('q'))
                
            Response(TestedPatchNum,:)= 0; %when subject press q to quit, response for current patch recorded as 0
            Confidence(TestedPatchNum,:) = 0;
            ResponseTime(TestedPatchNum, :) = 0;
            
                break
                
            end
       
        
    end %end trial loop
    
    % Now the wrap-up, we present the 'End of Run' text, save, and close the
    % screen
    Screen('FillRect', Exp.Cfg.win, Exp.Cfg.Color.gray);
    DrawFormattedText(Exp.Cfg.win,Exp.End_Sesh,'center','center',Exp.Cfg.Color.black);
    Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
    
    WaitSecs(wait_time_bw_screens); % Wait a amoment
    
    % ...And close the screen
    Screen('CloseAll');
    
    disp(['Thank you for your time, ' subj.initials '!']);
    Trialresponse = [subnumber CurrentTrialNumber QuestionNumber ImageFeature PatchSource PatchFeature Patch_location Response Confidence ResponseTime];
    save(save_path,'Trialresponse','-mat');
save(save_path, 'TrialResponse','-mat');
        disp('Data file saved!')    
    
    
catch ER
    ER.getReport
    sca;
end
