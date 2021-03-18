% NMRmulti2D: a GUI to plot one or more 2D NMR spectra. When plotting
% multiple 2D spectra, the F1 axis is shared between the spectra. Ideally
% only between 1-3 should be plotted side-by-side otherwise you'll have a
% very long figure. It is possible to magnify sections of any of the
% projections by using the provided button and clicking with the mouse to 
% define the line and the limits, in that order.

% Version 1.1 (01/28/2021)
% Recent changes:-
% Added functionality to plot magnified sections of the F1 axis/F2 axes

function NMRmulti2D
%% Create figure
h.f(1) = figure('units','normalized','position',[0.05,0.15,0.5,0.775],...
             'toolbar','none','menu','none');

%% Panels
h.panel(1) = uipanel('Title','F1','Units','normalized',...
                'FontSize',14,'FontWeight','bold',...
                'Position',[0.025,0.8,0.95,0.175]);
h.panel(2) = uipanel('Title','F2','Units','normalized',...
                'FontSize',14,'FontWeight','bold',...
                'Position',[0.025,0.45,0.45,0.35]);
h.panel(3) = uipanel('Title','2D Data','Units','normalized',...
                'FontSize',14,'FontWeight','bold',...
                'Position',[0.475,0.45,0.5,0.35]);

%% F1 Input
% Boxes
h.f1input(1) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.25,0.9,0.2,0.05],'String','Input F1 Folder','FontSize',14); 
h.f1input(2) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.5,0.9,0.1,0.05],'String','F1 Nucleus','FontSize',14); 

% Popup Menu
h.f1input(3) = uicontrol('Style','popupmenu','Units','normalized',...
                'Position',[0.5,0.875,0.1,0.05],'String',{'<HTML><sup>1</sup>H','<HTML><sup>13</sup>C','<HTML><sup>27</sup>Al','<HTML><sup>14</sup>N'},'Callback',@updatelimits); %F1 dropdown
% Buttons
h.f1input(4) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.9,0.2,0.05],'string','Browse Data','FontSize',14,...
                'callback',@p_browse_f1);
h.f1input(5) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.65,0.9,0.3,0.05],'string','Load F1','FontSize',14,...
                'callback',@p_load);    
     %% F1 Browse
    function p_browse_f1(varargin)

        f1name = uigetdir('~/','Choose the folder of the F1 dimension');
        if f1name == 0
            set(h.f1input(1),'String','Input 1D Folder')
        else
            set(h.f1input(1),'String',f1name)
        end

    end
    %% F1 Load
    function p_load(varargin)

        f1 = h.f1input(1).String;
        [x1, y1, ~] = brukimp1d(f1,1);
        DATA.x1 = x1;
        DATA.y1 = y1;
        disp('Loaded external projection in the F1 dimension')

        guidata(h.f(1),DATA);
    end
%% F2 Input
% Boxes
h.f2input(1) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.25,0.7,0.2,0.05],'String','Input F2 Dimension','FontSize',14);
h.f2input(2) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.2,0.6,0.1,0.05],'String','F2 Nucleus','FontSize',14); 
% Popup Menu
h.f2input(3) = uicontrol('Style','popupmenu','Units','normalized',...
                'Position',[0.2,0.575,0.1,0.05],'String',{'<HTML><sup>1</sup>H','<HTML><sup>13</sup>C','<HTML><sup>27</sup>Al','<HTML><sup>14</sup>N'},'Callback',@updatelimits); %F2 dropdown
set(h.f2input(3),'Value',3) % initial value of 27Al

% Buttons
h.f2input(4) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.7,0.2,0.05],'string','Browse Data','FontSize',14,...
                'callback',@p_browse_f2);
    %% F2 Browse
    function p_browse_f2(varargin)
        
%         [f2name, f2path] = uigetfile('/Users/leo/Dropbox (City College)/NMR Data/*.txt','Choose F2 Dimension');
%         f2name = strcat([f2path f2name]);
%         set(h.c(2),'String',f2name)
        
        f2name = uigetdir('~/','Choose the folder of the F2 dimension');
        if f2name == 0
            set(h.f2input(1),'String','Input 1D Folder')
        else
            set(h.f2input(1),'String',f2name)
        end

    end
%% 2D Input
h.input2D(1) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.725,0.7,0.2,0.05],'String','Input 2D Folder','FontSize',14); 
            
% Buttons
h.input2D(2) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.525,0.7,0.2,0.05],'string','Browse Data','FontSize',14,...
                'callback',@p_browsefolder);
    %% 2D Browse Callback
    function p_browsefolder(varargin)
       
        filefolder = uigetdir('~/','Choose the folder of the 2D experiment');
        if filefolder == 0
            set(h.input2D(1),'String','Input 2D Folder')
        else
            set(h.input2D(1),'String',filefolder)
        end
        
    end
%% Limits Panel
% Panel
% h.limits(1) = uipanel('Title','Limits','Units','normalized',...
%                 'FontSize',14,'FontWeight','bold',...
%                 'Position',[0.05,0.575,0.4,0.175]);
% Boxes
%F1
h.limits(1) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.175,0.85,0.1,0.035],'String','Lower','FontSize',14); %Lower String
h.limits(2) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.325,0.85,0.1,0.035],'String','Upper','FontSize',14); %Upper String
h.limits(3) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.175,0.825,0.1,0.035],'String','','FontSize',14); %Lower limit
h.limits(4) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.325,0.825,0.1,0.035],'String','','FontSize',14); %Upper limit
h.limits(5) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.0775,0.825,0.05,0.05],'String','F1 Limits','FontSize',14); %Label

%F2
h.limits(6) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.175,0.525,0.1,0.035],'String','Lower','FontSize',14); %Lower String
h.limits(7) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.325,0.525,0.1,0.035],'String','Upper','FontSize',14); %Upper String
h.limits(8) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.175,0.5,0.1,0.035],'String','','FontSize',14); %F2 Lower limit
h.limits(9) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.325,0.5,0.1,0.035],'String','','FontSize',14); %F2 Upper limit
h.limits(10) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.075,0.5,0.05,0.05],'String','F2 Limits','FontSize',14); %F2 Label

            
%% Populate Limits Boxes

    updatelimits(h.f1input(3));
    updatelimits(h.f2input(3));
     
    function [lowerlim, upperlim] = updatelimits(varargin)
        
        fninput = varargin{1};
        nucval = get(fninput, 'Value');

            if nucval == 1
                lowerlim = -2;
                upperlim = 18;
            elseif nucval == 2
                lowerlim = 0;
                upperlim = 200;
            elseif nucval == 3
                lowerlim = -40;
                upperlim = 120;
            elseif nucval == 4
                lowerlim = -200;
                upperlim = 200;
            end
            
            if fninput == h.f1input(3)
                set(h.limits(3),'String',num2str(lowerlim))
                set(h.limits(4),'String',num2str(upperlim))
            elseif fninput == h.f2input(3)
                set(h.limits(8),'String',num2str(lowerlim))
                set(h.limits(9),'String',num2str(upperlim))
            else
                disp('Problem with identifying nucleus inputs')
            end
    end

%% Threshold
% Boxes
h.T(1) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.55,0.52,0.4,0.05],'String','Threshold Factor','FontSize',14); %Threshold string
h.T(2) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.77,0.46,0.2,0.05],'String','% of maximum signal','FontSize',14); % % of maximum signal string
h.T(3) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.695,0.46,0.1,0.05],'String','','FontSize',14); %Slider edit box - threshold
% Slider
h.T(4) = uicontrol('Style','slider','Units','normalized',...
               'Position',[0.5,0.485,0.45,0.05],'Min',0.001,'Value',0.1,'SliderStep',[0.01001 0.1]);
% Live update of threshold text box
fun1 = @(~,e)set(h.T(3),'String',num2str(100*get(e.AffectedObject,'Value')));
addlistener(h.T(4), 'Value', 'PostSet',fun1);
set(h.T(3),'String',num2str(100*h.T(4).Value))
%% Number of Levels
% Boxes
h.levels(1) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.65,0.6125,0.2,0.05],'String','Number of Levels','FontSize',14); 
% Popup Menu
h.levels(2) = uicontrol('Style','popupmenu','Units','normalized',...
               'Position',[0.65,0.5825,0.2,0.05],'String',{'1','2','3','4','5','6','7','8','9','10'});
set(h.levels(2),'Value',6) % initial value of 6

        
%% Plot

h.p(1) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.35,0.45,0.075],'string','Plot First','FontSize',14,...
                'callback',@p_plot);
h.p(2) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.265,0.45,0.075],'string','Plot Next','FontSize',14,...
                'callback',@p_plotnext);
            
%% Magnification
h.mag(1) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.55,0.35,0.2,0.035],'String','Magnification Level','FontSize',14); % Magnification Box Text
h.mag(2) = uicontrol('Style','popupmenu','Units','normalized',...
                'Position',[0.75,0.335,0.15,0.05],'String',{'2x','4x','8x','16x','32x'},'Callback',@magnumber); % Magnification dropdown

            magnumber_initial = 2^h.mag(2).Value;
            
            DATA.magnumber = magnumber_initial;
            
h.mag(3) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.55,0.285,0.4,0.05],'string','Plot Magnified Section','FontSize',14,...
                'callback',@p_mag);    
            
%% Save
h.s(1) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.05,0.2,0.7,0.05],'String','Input Save File Name','FontSize',14); %Save name input
h.s(2) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.75,0.2,0.2,0.05],'string','Save','FontSize',14,...
                'callback',@p_save);
    %% Save Callback
    function p_save(varargin)
        savenmrfig(h.s(1).String)
                function savenmrfig(name)
                DATA = guidata(h.f(1));
                if exist('DATA.count')
                    count = DATA.count;
                else
                    count = 1;
                end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Here you can make the savelocation variable to suit 
                    % your computer

                        savelocation = '~';

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    filenamepdf=strcat(string(name), '.pdf');

                    orient(h.f(2),'landscape')
                    h.f(2).PaperSize = [3.5+6.15*(count),9];
                    h.f(2).Renderer='painters'; % ensures the figure is actually saved as a vector graphic
                    print(h.f(2),filenamepdf,'-dpdf','-fillpage')
                    movefile(filenamepdf,savelocation)
                    disp(strcat(["Figure saved to path: "], savelocation))
                end
    end
%% Logo and explanation            
h.exp(1) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.3,0.02,0.65,0.17],'String','Load F1 data first, using the browse function and the "Load F1" button. Then populated the F2 and 2D sections, tweaking the settings as desired. After desired settings have been set, push "Plot First". Then the F2 and 2D boxes can be repopulated with the next set of data and then push "Plot Next" to add it adjacent. This can be done ad nauseum. To magnify a section, use the button and click once to select the line, then twice more to define the limits. Code prepared by Leo Gordon, contact at lgordon@ccny.cuny.edu for questions.','FontSize',16);         
    try
        % Load in a background image and display it using the correct colors
        % The image used below, is in the Image Processing Toolbox.  If you do not have %access to this toolbox, you can use another image file instead.
        I=imread('./Messinger_Logo_Design.png');
        % This creates the 'background' axes
        ha = axes('units','normalized', ...
                    'position',[0.05 0.03 0.2 0.16]);
        % Move the background axes to the bottom
        uistack(ha,'bottom');

        hi = imagesc(I);
        % Turn the handlevisibility off so that we don't inadvertently plot into the axes again
        % Also, make the axes invisible
        set(ha,'handlevisibility','off', ...
                    'visible','off')  
    catch
    end
       
 %% Pushbutton callbacks

    function p_plot(varargin)
        
            DATA = guidata(h.f(1));
            
            fileinput = h.input2D(1).String;
            proc_no = 1;
            f2 = h.f2input(1).String;
            x1 = DATA.x1;
            y1 = DATA.y1;
            
            try
                if fileinput ~= "Input 2D Folder"
                    [Spectrum,~,xppm,yppm] = brukimp2d(fileinput, proc_no);
                else
                    disp('Please input folder')
                end
            catch
                disp('Error in reading file, please fix 2D folder input.')
                h.input2D(1).String = 'Input 2D folder';
                return
            end
            
            
            if f2 ~= "Input F2 Dimension"
                [x2, y2, ~] = brukimp1d(f2,1);
                disp('Loaded external projection in the F2 dimension')
            else
                disp('Please input filepath for the F2 dimension')
            end

        factor = h.T(4).Value;
        xc = xppm;
        yc = yppm;  

        logicalll = Spectrum>0; % Positive values only
        thres = logicalll.*Spectrum;
        
        thresmax = max(thres(:));      
        thresmin = factor*thresmax;
        clevels = h.levels(2).Value;
        thresvec=nan(1,clevels);
        for i=1:clevels
            thresvec(i)=thresmin*((thresmax/thresmin)^(1/clevels))^i;
        end
        
        f2l = h.limits(8).String;
        f2u = h.limits(9).String;
        f1l = h.limits(3).String;
        f1u = h.limits(4).String;
        
        F2_nuc_num = h.f2input(3).Value;
        F1_nuc_num = h.f1input(3).Value;
        
        if (~isempty(f2u) && ~isempty(f2l))
            F2upper = str2double(f2u);
            F2lower = str2double(f2l);
            [F2ax, xstr, F2ticks] = nucleus_choice(F2_nuc_num,F2lower,F2upper);
        else
            disp('Incomplete F2 axis limits, using defaults')
            [F2ax, xstr, F2ticks] = nucleus_choice(F2_nuc_num);
        end
        
        if (~isempty(f1u) && ~isempty(f1l))
            F1upper = str2double(f1u);
            F1lower = str2double(f1l);
            [F1ax, ystr, F1ticks] = nucleus_choice(F1_nuc_num,F1lower,F1upper);
        else
            disp('Incomplete F1 axis limits, using defaults')
            [F1ax, ystr, F1ticks] = nucleus_choice(F1_nuc_num);            
        end

        h.f(2) = figure('Units','inches','Position',[0 0 10 10]);
        conpos = [2.5 1.5 6 6];
            
        conposy = conpos;
        conposy(1) = conpos(1)-2.075;
        conposy(3) = 2;
          
        conposx = conpos;
        conposx(2) = conpos(2)+6.075;
        conposx(4) = 2;
        
        a1 = axes(h.f(2),'Units','inches','Position',conpos,'box','on','xtick',[],'ytick',[],'LineWidth',2);  %Under axis, used as background color and box;        
        axc = axes(h.f(2),'Units','inches','Position',conpos);
        propLink = linkprop([a1 axc],'position');      
        contour(xc,yc,thres,thresvec(:),'.-k','LineWidth',1.5);
        axis([F2ax F1ax])
        axc.Box='off';
        axc.Color='none';
        axc.XDir='reverse';
        axc.YDir='reverse';
        axc.YAxisLocation='right';
        xlabel(xstr)
        ylabel(ystr)
        axc.FontWeight = 'bold';
        axc.LineWidth = 2;
        axc.FontSize = 18;
        axc.XAxis.TickDirection='out';
        axc.YAxis.TickDirection='out';
        xticks(F2ticks)
        yticks(F1ticks)
        axc.XMinorTick = 'on';
        axc.YMinorTick = 'on';

        % [left bottom width height]
        % x -> [leftsame bottom+height widthsame height=0.2]
        % y -> [left-width bottomsame width=0.2 heightsame]

        ax1=axes(h.f(2),'Units','inches','Position',conposy);
        H1 = plot(ax1,y1,x1,'-k','LineWidth',3);
        ax1.XDir='reverse';
        ax1.YDir='reverse';
        ax1.YLim=F1ax;
        ax1.Box='off';
        ax1.XAxis.Visible='off';
        ax1.YAxis.Visible='off';

        ax2=axes(h.f(2),'Units','inches','Position',conposx);
        H2 = plot(ax2,x2,y2,'-k','LineWidth',3);
        ax2.XDir='reverse';
        ax2.XLim=F2ax;
        ax2.Box='off';
        ax2.XAxis.Visible='off';
        ax2.YAxis.Visible='off';
        
        DATA.axc = axc;
        DATA.F1ax = F1ax;
        DATA.F2ax = F2ax;
        DATA.figure = h.f(2);        
        DATA.ystr = ystr;
        DATA.F1ticks = F1ticks;
        DATA.H1 = H1;
        DATA.H2 = H2;
        DATA.ax1 = ax1;
        DATA.ax2 = ax2;
        guidata(h.f(1),DATA);

        disp('Done Plotting')

    end

    function p_plotnext(varargin)
        DATA = guidata(h.f(1));
        
        datafig = DATA.figure;
        F1ax = DATA.F1ax;
        F2axold = DATA.F2ax;
        axc = DATA.axc;
        ystr = DATA.ystr;
        F1ticks = DATA.F1ticks;
            
            fileinput = h.input2D(1).String;
            proc_no = 1;
            f2 = h.f2input(1).String;
            
            try
                if fileinput ~= "Input 2D Folder"
                    [Spectrum,~,xppm,yppm] = brukimp2d(fileinput, proc_no);
                else
                    disp('Please input folder')
                end
            catch
                disp('Error in reading file, please fix 2D folder input.')
                h.input2D(1).String = 'Input 2D folder';
                return
            end
            
            
            if f2 ~= "Input F2 Dimension"
                [x2, y2, ~] = brukimp1d(f2,1);
                disp('Loaded external projection in the F2 dimension')
            else
                disp('Please input filepath for the F2 dimension')
            end

        factor = h.T(4).Value;
        xc = xppm;
        yc = yppm;  

        logicalll = Spectrum>0; % Positive values only
        thres = logicalll.*Spectrum;
        
        thresmax = max(thres(:));      
        thresmin = factor*thresmax;
        clevels = h.levels(2).Value;
        thresvec=nan(1,clevels);
        for i=1:clevels
            thresvec(i)=thresmin*((thresmax/thresmin)^(1/clevels))^i;
        end
        
        f2l = h.limits(8).String;
        f2u = h.limits(9).String;
        nuc_num = h.f2input(3).Value;
        
        if (~isempty(f2u) && ~isempty(f2l))
            F2upper = str2double(f2u);
            F2lower = str2double(f2l);
            [F2ax, xstr, F2ticks] = nucleus_choice(nuc_num,F2lower,F2upper);
        else
            disp('Incomplete F2 axis limits, using defaults')
            [F2ax, xstr, F2ticks] = nucleus_choice(nuc_num);
        end

        n = length(findall(datafig,'type','axes'));
        n1 = length(findall(datafig,'type','line'));
        i=n1-1;
        if n == 4 % find number of axes instead
                count = 1;
            axc.YAxis.TickLabels = [];
            axc.XTick(axc.XTick==F2axold(1))=[];
        else
            count = DATA.count;
            F2axold2 = DATA.F2ax2;
            axc2 = DATA.axc2;
            axc2.YAxis.TickLabels = [];
            axc2.XTick(axc2.XTick==F2axold2(1))=[];        
        end
        
        conpos = [2.5+(6.15*(count)) 1.5 6 6];  
        conposx = conpos;
        conposx(2) = conpos(2)+6.075;
        conposx(4) = 2;
%look into turning off right ticks after each
% set a1 visibility to off? findall vs findobj
        a1 = axes(h.f(2),'Units','inches','Position',conpos,'box','on','xtick',[],'ytick',[],'LineWidth',2);  %Under axis, used as background color and box;        
        
        axc2 = axes(datafig,'Units','inches','Position',conpos);  
        propLink = linkprop([a1 axc2],'position');      
        contour(axc2,xc,yc,thres,thresvec(:),'.-k','LineWidth',1.5);        
        xlabel(axc2,xstr)
        ylabel(axc2,ystr)
        axis(axc2,[F2ax F1ax])
        axc2.Box='off';
        axc2.Color='none';
        axc2.XDir='reverse';
        axc2.YDir='reverse';
        axc2.YAxisLocation='right';

        axc2.FontWeight = 'bold';
        axc2.LineWidth = 2;
        axc2.FontSize = 18;
        axc2.XAxis.TickDirection='out';
        axc2.YAxis.TickDirection='out';
        xticks(axc2,F2ticks)
        yticks(axc2,F1ticks)
        axc2.XMinorTick = 'on';
        axc2.YMinorTick = 'on';

        % [left bottom width height]
        % x -> [leftsame bottom+height widthsame height=0.2]
        % y -> [left-width bottomsame width=0.2 heightsame]

        ax22=axes(datafig,'Units','inches','Position',conposx);
        H(i) = plot(ax22,x2,y2,'-k','LineWidth',3);
        ax22.XDir='reverse';
        ax22.XLim=F2ax;
        ax22.Box='off';
        ax22.XAxis.Visible='off';
        ax22.YAxis.Visible='off';
        
        count = count+1;
        datafig.Position = [0 0 2.5+7*(count) 10];
        
        DATA.axc2 = axc2;
        DATA.F2ax2 = F2ax;
        DATA.n = n;
        DATA.count = count;
        DATA.H(i) = H(i);
        DATA.ax22(i) = ax22;
        guidata(h.f(1),DATA);
        
        disp('Done Plotting')

    end

    function magnumber(varargin)
            magupdate = h.mag(2).Value;

            magnumber = 2^magupdate;
            
            DATA.magnumber = magnumber;
            guidata(h.f(1),DATA);
    end

    function DATA = p_mag(varargin)
        
        DATA = guidata(h.f(1));
        H1 = DATA.H1;
        H2 = DATA.H2;
        try 
            H =  DATA.H;
            set(H, 'ButtonDownFcn', {@LineSelected, H})
        catch
        end
        
        try
            set(H2, 'ButtonDownFcn', {@LineSelected, H2})
        catch
        end
        
        try
            set(H1, 'ButtonDownFcn', {@LineSelected1, H1})
        catch
        end
        return
        
        function LineSelected(ObjectH, EventData, Hfn)
            
            clear ax
            DATA = guidata(h.f(1));

            if Hfn == H2
                ax = DATA.ax2;
            elseif Hfn == H
                try
                    ancest = ancestor(gcbo,'axes');
                    ax = ancest; % need to find this index somehow
                catch
                end
            end
            magnification = DATA.magnumber;

            [magx, ~] = ginput_ax(2);
            xdata = get(ObjectH, 'XData')';
            ydata = get(ObjectH, 'YData')';
            colour = get(ObjectH, 'Color');
            
            if magx(2)>magx(1)
                istart = find(xdata>=magx(2),1,'last'); %gets higher ppm value index
                iend = find(xdata>=magx(1),1,'last'); % gets lower ppm value index
            else
                istart = find(xdata>=magx(1),1,'last'); %gets higher ppm value index
                iend = find(xdata>=magx(2),1,'last'); % gets lower ppm value index
            end
            
            baseline = ydata(istart);
            
            xmag = xdata(istart:iend);
            ymag = ydata(istart:iend).*magnification;
            
            magbaseline = min(ymag);
%             ymag = ymag-min(ymag);
%             ymag = ymag./max(ymag);
            ymag = (ymag)+baseline-magbaseline+0.1;
            hold on
            plot(ax,xmag,ymag,'-','Color',colour,'LineWidth',2)
            
%             label?
            text(xmag(end)-0.5,ymag(end),strcat(['\times',num2str(magnification)]),'FontSize',12,'FontWeight','bold')
            hold off
            return
        end
        
        function LineSelected1(ObjectH, EventData, H)
            
            DATA = guidata(h.f(1));
            ax = DATA.ax1;
            magnification = DATA.magnumber;

            [~, magy] = ginput_ax(2);
            xdata = get(ObjectH, 'XData')';
            ydata = get(ObjectH, 'YData')';
            colour = get(ObjectH, 'Color');
            
            if magy(2)>magy(1)
                istart = find(ydata>=magy(2),1,'last'); %gets higher ppm value index
                iend = find(ydata>=magy(1),1,'last'); % gets lower ppm value index
            else
                istart = find(ydata>=magy(1),1,'last'); %gets higher ppm value index
                iend = find(ydata>=magy(2),1,'last'); % gets lower ppm value index
            end
            
            baseline = xdata(istart);
            
            xmag = xdata(istart:iend).*magnification;
            ymag = ydata(istart:iend);
            
            magbaseline = min(xmag);
% %             ymag = ymag-min(ymag);
% %             ymag = ymag./max(ymag);
            xmag = (xmag)+baseline-magbaseline+0.1;
            hold on
            plot(ax,xmag,ymag,'-','Color',colour,'LineWidth',2)
            
%             label?
            text(xmag(end),ymag(end)-0.5,strcat(['\times',num2str(magnification)]),'FontSize',12,'FontWeight','bold','Rotation',90)
            hold off
            return
        end
        % H and H1 now can assign to lines by clicking, use these inputs to
        % do the magnification
    end
%% Called Functions

%%%%%%
    function [Spectrum, Params,xppm,yppm] = brukimp2d(fileinput, proc_no)

                input_file = strcat(fileinput,'/pdata/',string(proc_no),'/2rr');
                Acqus_2D_H_File_Path = string([fileinput '/acqus']);
                Procs_2D_H_File_Path = strcat(fileinput, '/pdata/', string(proc_no), '/procs');
                Acqus_2D_C_File_Path = string([fileinput '/acqu2s']);
                Procs_2D_C_File_Path = strcat(fileinput, '/pdata/', string(proc_no), '/proc2s');

                fid = fopen(input_file, 'rb');
                if fid < 1
                    error('File not found %s\n', input_file);
                else
                    Spectrum_2D = fread(fid, 'int');
                end
                fclose(fid);

                fid_aqus = fopen(Acqus_2D_H_File_Path, 'r');
                fid_procs = fopen(Procs_2D_H_File_Path, 'r');
                if fid_aqus < 1 || fid_procs < 1
                    fclose(fid_aqus);
                    fclose(fid_procs);
                    error('Could not open %s or %s\n', Acqus_2D_H_File_Path, Procs_2D_H_File_Path);
                else
                    [H_OBS, H_CAR, H_Error_aqus] = Get_Bruker_Info_1D_Acqus(fid_aqus);
                    [H_SF, H_SW, H_Length, H_Error_proc] = Get_Bruker_Info_1D_Procs(fid_procs);
                    if ~isempty(H_Error_aqus) || ~isempty(H_Error_proc)
                        fclose(fid_aqus);
                        fclose(fid_procs);
                        error('Something went wrong with the params in %s or %s\n', Acqus_2D_H_File_Path, Procs_2D_H_File_Path);
                    end
                end
                fclose(fid_aqus);
                fclose(fid_procs);

                Params.xOBS = H_OBS;
                Params.xCAR = H_CAR;
                Params.xSW = H_SW;
                Params.xSF = H_SF;

                fid_acqus = fopen(Acqus_2D_C_File_Path, 'r');
                fid_procs = fopen(Procs_2D_C_File_Path, 'r');
                if fid_acqus < 1 || fid_procs < 1
                    fclose(fid_acqus);
                    fclose(fid_procs);
                    error('Could not open  %s or %s \n', Acqus_2D_C_File_Path, Procs_2D_C_File_Path);
                else
                    [C_OBS, C_CAR, C_Error_aqus] = Get_Bruker_Info_1D_Acqus(fid_acqus);
                    [C_SF, C_SW, C_Length, C_Error_proc] = Get_Bruker_Info_1D_Procs(fid_procs);
                    if ~isempty(C_Error_aqus) || ~isempty(C_Error_proc)
                        fclose(fid_acqus);
                        fclose(fid_procs);
                        error('Something went wrong with the params in %s or %s\n', Acqus_2D_C_File_Path, Procs_2D_C_File_Path);
                    end
                end
                fclose(fid_acqus);
                fclose(fid_procs);

                Params.yOBS = C_OBS;
                Params.yCAR = C_CAR;
                Params.ySW = C_SW;
                Params.ySF = C_SF;

                Spectrum = reshape(Spectrum_2D, H_Length, C_Length);
                Spectrum = Spectrum';

                xaxcen=Params.xCAR*Params.xOBS-((Params.xSF-Params.xOBS)*1000000); % why don't I just read in SFO1? not same?
                xaxmin=xaxcen-Params.xSW/2;
                xaxmax=xaxcen+Params.xSW/2;
                xaxlen=(xaxmax-xaxmin)/(length(Spectrum)-1);
                xaxhz=(xaxmin:xaxlen:xaxmax);
                xppm=xaxhz/Params.xOBS;
                xppm=sort(xppm,'descend');

                yaxcen=Params.yCAR*Params.yOBS-((Params.ySF-Params.yOBS)*1000000); % why don't I just read in SFO1?
                yaxmin=yaxcen-Params.ySW/2;
                yaxmax=yaxcen+Params.ySW/2;
                yaxlen=(yaxmax-yaxmin)/(size(Spectrum,1)-1);
                yaxhz=(yaxmin:yaxlen:yaxmax);
                yppm=yaxhz/Params.yOBS;
                yppm=sort(yppm,'descend');
            end
%%%%%%

%%%%%%
    function [xppm, Spectrum, Params] = brukimp1d(fileinput, proc_no)

        input_file = strcat(fileinput,'/pdata/',string(proc_no),'/1r');
        Acqus_1D_H_File_Path = string([fileinput '/acqus']);
        Procs_1D_H_File_Path = strcat(fileinput, '/pdata/', string(proc_no), '/procs');

        fid = fopen(input_file, 'rb');
        if fid < 1
            error('File not found %s\n', input_file);
        else
            Spectrum_1D = fread(fid, 'int');
        end
        fclose(fid);

        fid_aqus = fopen(Acqus_1D_H_File_Path, 'r');
        fid_procs = fopen(Procs_1D_H_File_Path, 'r');
        if fid_aqus < 1 || fid_procs < 1
            fclose(fid_aqus);
            fclose(fid_procs);
            error('Could not open %s or %s\n', Acqus_1D_H_File_Path, Procs_1D_H_File_Path);
        else
            [H_OBS, H_CAR, H_Error_aqus] = Get_Bruker_Info_1D_Acqus(fid_aqus);
            [H_SF, H_SW, H_Length, H_Error_proc] = Get_Bruker_Info_1D_Procs(fid_procs);
            if ~isempty(H_Error_aqus) || ~isempty(H_Error_proc)
                fclose(fid_aqus);
                fclose(fid_procs);
                error('Something went wrong with the params in %s or %s\n', Acqus_1D_H_File_Path, Procs_1D_H_File_Path);
            end
        end
        fclose(fid_aqus);
        fclose(fid_procs);

        Params.xOBS = H_OBS;
        Params.xCAR = H_CAR;
        Params.xSW = H_SW;
        Params.xSF = H_SF;

        Spectrum = Spectrum_1D;
        Spectrum = Spectrum-min(Spectrum);
        Spectrum = Spectrum./max(Spectrum);

        xaxcen=Params.xCAR*Params.xOBS-((Params.xSF-Params.xOBS)*1000000); % why don't I just read in SFO1? not same?
        xaxmin=xaxcen-Params.xSW/2;
        xaxmax=xaxcen+Params.xSW/2;
        xaxlen=(xaxmax-xaxmin)/(H_Length-1);
        xaxhz=(xaxmin:xaxlen:xaxmax);
        xppm=xaxhz/Params.xOBS;
        xppm=sort(xppm,'descend');
        
    end
%%%%%%

%%%%%%
    function [SF, SW, Length, Error] = Get_Bruker_Info_1D_Procs(fid)

            SW = 0;
            Length = 0;
            SF = 0;

            tline = fgetl(fid);
            Satisfied = false;
            while ~Satisfied
                    if ~isempty(strfind(tline, '##$SW_p= '))
                        tline = strrep(tline, '##$SW_p= ', '');
                        SW = str2double(tline);
                    end
                    if ~isempty(strfind(tline, '##$SI= '))
                        tline = strrep(tline, '##$SI= ', '');        
                        Length = str2double(tline);
                    end
                    if ~isempty(strfind(tline, '##$SF= '))
                        tline = strrep(tline, '##$SF= ', '');
                        SF = str2double(tline);
                    end
                tline = fgetl(fid);
                    if ~ischar(tline) || (SW~=0 && Length~= 0 && SF~=0)
                        Satisfied = true;
                    end
            end


                if (SW~=0 && Length ~= 0)
                    Error = '';
                else
                    Error = 'Could not find all the parameters from the aqcus file';
                end
        end
%%%%%%

%%%%%%
    function [OBS, CAR, Error] = Get_Bruker_Info_1D_Acqus(fid)

            OBS = nan;
            CAR = nan;
            O1 = nan;

            tline = fgetl(fid);
            Satisfied = false;
            while ~Satisfied
                if ~isempty(strfind(tline, '##$O1= '))
                    tline = strrep(tline, '##$O1= ', '');
                    O1 = str2double(tline);
                end
                if ~isempty(strfind(tline, '##$BF1= '))
                    tline = strrep(tline, '##$BF1= ', '');
                    OBS = str2double(tline);
                end
                tline = fgetl(fid);
                if ~isnan(OBS) && ~isnan(O1)
                    Satisfied = true;
                end
            end

            if (OBS~=0 && O1~=0)
                CAR = O1/OBS;
                Error = '';
            elseif O1==0
                CAR = 1e-30;
                Error = '';
            else
                Error = 'Could not find all the parameters from the aqcus file';
            end
        end
%%%%%%
    function [axlims,labelstr,axticks] = nucleus_choice(varargin)
        nuc_num = varargin{1};
        if nargin == 3
            F2lower = varargin{2};
            F2upper = varargin{3};
        end
        if nuc_num == 1
%                 F2_nuc = 1H;
            try
                axlims = [F2lower F2upper];
            catch
                axlims = [-2 18];
            end
            axticks = [-100:2:300];
            labelstr = '^{1}H chemical shift (ppm)';
        elseif nuc_num == 2
%                 F2_nuc = 13C;
            try
                axlims = [F2lower F2upper];
            catch
                axlims = [0 200];
            end
            axticks = [-100:20:300];
            labelstr = '^{13}C chemical shift (ppm)';
        elseif nuc_num == 3
%                 F2_nuc = 27Al;
            try
                axlims = [F2lower F2upper];
            catch
                axlims = [-40 120];
            end
            axticks = [-100:20:300];
            labelstr = '^{27}Al shift (ppm)';
        elseif nuc_num == 4
%                 F2_nuc = 14N;
            try
                axlims = [F2lower F2upper];
            catch
                axlims = [-80 80];
            end
            axticks = [-100:20:300];
            labelstr = '^{14}N shift (ppm)';
        end
    end
%%%%%%
    function varargout = ginput_ax(n)

    k = 0;
    xy = zeros(n,2);
    hf = h.f(2);    
%     hf = get(ha,'parent');
    ha = hf.CurrentAxes;
    figure(hf);
    set(hf,'WindowButtonMotionFcn',@changepointer)
    set(ha,'ButtonDownFcn',@getpoints)
    hp = get(ha,'children');
    ht = get(hp,'hittest');
    set(hp,'hittest','off')
    axlim = get(ha,'Position');
    fglim = get(hf,'Position');
    x1 = axlim(1)*fglim(3) + fglim(1);
    x2 = (axlim(1)+axlim(3))*fglim(3) + fglim(1);
    y1 = axlim(2)*fglim(4) + fglim(2);
    y2 = (axlim(2)+axlim(4))*fglim(4) + fglim(2);
    waitfor(hf,'WindowButtonMotionFcn',[])
    if iscell(ht)
        for jj=1:length(ht)
            set(hp(jj),'hittest',ht{jj})
        end
    else
        set(hp,'hittest',ht)
    end
    if nargout==2
        varargout{1} = xy(:,1);
        varargout{2} = xy(:,2);
    else
        varargout{1} = xy;
    end
          function changepointer(~,~)
              pntr = get(0,'PointerLocation');
              if pntr(1)>x1 && pntr(1)<x2 && pntr(2)>y1 && pntr(2)<y2
                  set(hf,'Pointer','crosshair')
              else
                  set(hf,'Pointer','arrow')
              end
          end
          function getpoints(hObj,~,~)
              cp = get(hObj,'CurrentPoint');
              k = k+1;
              xy(k,:) = cp(1,1:2);
              if k==n
                  set(hf,'Pointer','arrow')
                  set(hf,'WindowButtonMotionFcn',[])
                  set(ha,'ButtonDownFcn',[])
              end
          end
          return
    end
%%% Function End %%%
end

