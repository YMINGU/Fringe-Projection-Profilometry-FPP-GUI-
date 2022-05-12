classdef FPP_YG_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        TabGroup                      matlab.ui.container.TabGroup
        CalibrationTab                matlab.ui.container.Tab
        KxyEditField                  matlab.ui.control.EditField
        KxyLabel                      matlab.ui.control.Label
        TextArea                      matlab.ui.control.TextArea
        TextAreaLabel                 matlab.ui.control.Label
        LoadButton                    matlab.ui.control.Button
        CalibrationButton             matlab.ui.control.Button
        UIAxes2                       matlab.ui.control.UIAxes
        TopographyCalculationTab      matlab.ui.container.Tab
        InstructionsTextArea          matlab.ui.control.TextArea
        InstructionsTextAreaLabel     matlab.ui.control.Label
        StartButton                   matlab.ui.control.Button
        BrowseReferencePatternButton  matlab.ui.control.Button
        BrowseFringePatternButton     matlab.ui.control.Button
        UIAxes                        matlab.ui.control.UIAxes
    end

    
    methods (Access = private)
        
        function M = p_shift(~,path)
            fp = fullfile(path, '*.jpg'); 
            img = dir(fp);
            Y = 0;
            X = 0;
            n = 3;
            
            for k = 1:3
                fpName = img(k).name;
                imgName = fullfile(path, fpName);
                II = imread(imgName);
                if size(II,3)==3
                    II = rgb2gray(II);
                end 
                II=double(II);
                sigma = 2*(k-1)*pi/n;
                Y = Y + sin(sigma) * II;
                X = X + cos(sigma) * II;
            end
            M = (atan2(-Y,X));
            
        end
        
        function [Iunwrapped] = itoh_unwrap_c(~,Iwrapped)
            [m,n] = size(Iwrapped);
            Iunwrapped = Iwrapped;
            
            for i=1:n
                Iunwrapped(end:-1:1,i) = unwrap(Iunwrapped(end:-1:1,i));
            end
            
            for i=1:m
                Iunwrapped(i,:) = unwrap(Iunwrapped(i,:));
            end
            
        end
        
        function uphi_ref = unfold(~,phi)
            [w,h] = size(phi);
            uphi_ref = phi;
            for i=1:h
                uphi_ref(:,i) = unwrap(uphi_ref(:,i));
            end
            for j=1:w
                uphi_ref(j,:) = unwrap(uphi_ref(j,:));
            end
            
            
        end
        
        
        function P = fitting(~,M)
            [w,h] = size(M);
            x = 1:h;
            x2 = (1:w)';
            
            parfor i = 1:w
                A = polyfit(x,M(i,:),2);
                y = A(1).*x.^2+A(2).*x.^1+A(3);
                M(i,:) = y;
            end
            
            parfor j = 1:h
                A = polyfit(x2,M(:,j),2);
                y = A(1).*x2.^2+A(2).*x2.^1+A(3);
                M(:,j) = y;
        
        
            end
            P=M;
        end
        
            function phi_u=unwrap_ref(app,phi,phi_ref)
            phi = phi +2*pi*round((phi_ref-phi)/(2*pi));
            phi2 =fitting(app,phi);
            dphi = phi-phi2;
            dphi = F_filter(app,dphi);
            uphi = dphi + phi2;
            phi_u = uphi;
            
            
        end
        
        
        function uphi=F_filter(~,phi)
            [w,h]=size(phi);
            fyy=32;
            fxx=5;
            rx=20;
            ry=5;
            cY=w/2+1;
            cX=h/2+1;
            FT=fftshift(fft2(phi));
            
            for k = 1:3
                fx = k*fxx; 
                fy = k*fyy;
                
                FT(floor(cY-ry):ceil(cY+ry),floor(cX-fx-rx):ceil(cX-fx+rx))=0;
                FT(floor(cY+fy-ry):ceil(cY+fy+ry),floor(cX-rx):ceil(cX+rx))=0;
                FT(floor(cY-ry):ceil(cY+ry),floor(cX+fx-rx):ceil(cX+fx+rx))=0;
                FT(floor(cY+fy-ry):ceil(cY+fy+ry),floor(cX-rx):ceil(cX+rx))=0;
            end

            FT(1:ry,floor(cX-rx):ceil(cX+rx))=0;
            FT(end-ry:end,floor(cX-rx):ceil(cX+rx))=0;
            FT(floor(cY-ry):ceil(cY+ry), 1:rx)=0;
            FT(floor(cY-ry):ceil(cY+ry), end-rx:end)=0;
            
            uphi = real(ifft2(ifftshift(FT)));
       end

            
            
            
            
        end
        
        
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            clear;
            clc;
            format long;
        end

        % Button pushed function: BrowseFringePatternButton
        function BrowseFringePatternButtonPushed(app, event)

            global filename
            [~,Pnameh]=uigetfile('*.jpg');
            filename=Pnameh;
        end

        % Button pushed function: BrowseReferencePatternButton
        function BrowseReferencePatternButtonPushed(app, event)
            global filenameR
            [~,PnamehR]=uigetfile('*.jpg');
            filenameR=PnamehR;
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            for i=1:3
                global filename
                global filenameR
                global Cal_k;
                path = filename;
                path2 = filenameR;
                phi_t = p_shift(app,path);
                phi_r = p_shift(app,path2);
                uphi_rr = unfold(app,phi_r);
                uphi_r=fitting(app,uphi_rr);
                uphi_t = unwrap_ref(app,phi_t,uphi_r);
                delta=cell(3,1);
                delta{i,1} =F_filter(app,uphi_t-uphi_r);
                z=cell(3,1);
                z{i,1} = cell2mat(delta(i,1)).*Cal_k;
                height_c_rxy = cell2mat(z(i,1));
                aheight=(1/2)*abs(height_c_rxy);
                mesh(app.UIAxes,aheight)
                view(app.UIAxes,[-144,88]);
                
            end
        end

        % Button pushed function: CalibrationButton
        function CalibrationButtonPushed(app, event)
            global filenameL
            n=10;
            H=1000;
            phi = cell(n,1);
            uphi = cell(n,1);
            for i = 1:n
                path=[filenameL,num2str(i),''];
                phi_i= p_shift(app,path);
                phi{i,1} = phi_i;
            end
            uphi_r = itoh_unwrap_c(app,cell2mat(phi(n,1)));
            uphi_r=fitting(app,uphi_r);
            [w,h] = size(uphi_r);
            parfor i = 1:n
            U = unwrap_ref(app,cell2mat(phi(i,1)),uphi_r);
            uphi{i,1} = U;  
            end
            
            Kk = zeros(w,h,2);
            y = zeros(n,1);
            yy = zeros(n,1);
            
            for i = 1
                for j = 1440
                    for z = 1:n
                        M = cell2mat(uphi(z,1));
                        y(z,1) = M(i,j);
                    end
                    x = linspace(0,n-1,n)';
                    K = polyfit(x,y,1);
                    Kk(i,j,1) = K(1);
                    Kk(i,j,2) = K(2);
                end 
            end
            scatter(app.UIAxes2,y,x)
            K=polyfit(y,x,1);
            
            global Cal_k
            Cal_k=K(1);
            app.KxyEditField.Value=num2str(Cal_k);
            
            
            
            
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            global filenameL
            [~,Cnameh]=uigetfile('*.jpg');
            filenameL=Cnameh;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 795 650];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 795 650];

            % Create CalibrationTab
            app.CalibrationTab = uitab(app.TabGroup);
            app.CalibrationTab.Title = 'Calibration';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.CalibrationTab);
            title(app.UIAxes2, 'Calibration of K(x,y) in the FPP System')
            xlabel(app.UIAxes2, 'Location, in mm')
            ylabel(app.UIAxes2, 'Phase, in rad')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Position = [51 177 688 428];

            % Create CalibrationButton
            app.CalibrationButton = uibutton(app.CalibrationTab, 'push');
            app.CalibrationButton.ButtonPushedFcn = createCallbackFcn(app, @CalibrationButtonPushed, true);
            app.CalibrationButton.Position = [51 97 165 22];
            app.CalibrationButton.Text = 'Calibration';

            % Create LoadButton
            app.LoadButton = uibutton(app.CalibrationTab, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Position = [51 128 165 22];
            app.LoadButton.Text = 'Load';

            % Create TextAreaLabel
            app.TextAreaLabel = uilabel(app.CalibrationTab);
            app.TextAreaLabel.HorizontalAlignment = 'right';
            app.TextAreaLabel.Position = [297 93 55 22];
            app.TextAreaLabel.Text = 'Text Area';

            % Create TextArea
            app.TextArea = uitextarea(app.CalibrationTab);
            app.TextArea.Editable = 'off';
            app.TextArea.Position = [367 57 372 93];
            app.TextArea.Value = {'1. Click "load" to load the folder path of the calibration photo'; '    (Click "0.jpg" to locate the folder)'; '2. Click "Calibration" to start Calibration'; '3. Calibrated value would be shown at the K(x,y) tab. '; 'FYI: Remeber to first do calibration before topography calucation.'; 'Copyrights by Yuming Gu (yug52@pitt.edu)'};

            % Create KxyLabel
            app.KxyLabel = uilabel(app.CalibrationTab);
            app.KxyLabel.HorizontalAlignment = 'right';
            app.KxyLabel.Position = [51 67 35 22];
            app.KxyLabel.Text = 'K(x,y)';

            % Create KxyEditField
            app.KxyEditField = uieditfield(app.CalibrationTab, 'text');
            app.KxyEditField.Position = [101 67 115 22];

            % Create TopographyCalculationTab
            app.TopographyCalculationTab = uitab(app.TabGroup);
            app.TopographyCalculationTab.Title = 'Topography Calculation';

            % Create UIAxes
            app.UIAxes = uiaxes(app.TopographyCalculationTab);
            title(app.UIAxes, 'Height Topography of 3D Pattern')
            xlabel(app.UIAxes, 'Horizontal Pixel')
            ylabel(app.UIAxes, 'Vertical Pixel')
            zlabel(app.UIAxes, '                Height, in cm')
            app.UIAxes.FontSize = 14;
            app.UIAxes.Position = [43 149 681 465];

            % Create BrowseFringePatternButton
            app.BrowseFringePatternButton = uibutton(app.TopographyCalculationTab, 'push');
            app.BrowseFringePatternButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseFringePatternButtonPushed, true);
            app.BrowseFringePatternButton.Tag = 'OpenFile';
            app.BrowseFringePatternButton.Position = [51 118 165 22];
            app.BrowseFringePatternButton.Text = 'Browse Fringe Pattern';

            % Create BrowseReferencePatternButton
            app.BrowseReferencePatternButton = uibutton(app.TopographyCalculationTab, 'push');
            app.BrowseReferencePatternButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseReferencePatternButtonPushed, true);
            app.BrowseReferencePatternButton.Tag = 'OpenFile';
            app.BrowseReferencePatternButton.Position = [51 82 165 22];
            app.BrowseReferencePatternButton.Text = 'Browse Reference Pattern';

            % Create StartButton
            app.StartButton = uibutton(app.TopographyCalculationTab, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.Position = [51 46 165 22];
            app.StartButton.Text = 'Start';

            % Create InstructionsTextAreaLabel
            app.InstructionsTextAreaLabel = uilabel(app.TopographyCalculationTab);
            app.InstructionsTextAreaLabel.HorizontalAlignment = 'right';
            app.InstructionsTextAreaLabel.Position = [263 82 68 22];
            app.InstructionsTextAreaLabel.Text = 'Instructions';

            % Create InstructionsTextArea
            app.InstructionsTextArea = uitextarea(app.TopographyCalculationTab);
            app.InstructionsTextArea.Editable = 'off';
            app.InstructionsTextArea.Position = [338 46 425 94];
            app.InstructionsTextArea.Value = {'1. Click "Browse Fringe Pattern" to load the folder of fringe pattern.'; '2. Click "Browse Reference Pattern" to load the folder of Reference Pattern.'; '3. Click"Start" to discover what will happened!'; ''; 'Copyrights by Yuming Gu (yug52@pitt.edu)'};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = FPP_YG_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end