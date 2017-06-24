function varargout = JTK(varargin)
% JTK MATLAB code for JTK.fig
%      JTK, by itself, creates a new JTK or raises the existing
%      singleton*.
%
%      H = JTK returns the handle to a new JTK or the handle to
%      the existing singleton*.
%
%      JTK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JTK.M with the given input arguments.
%
%      JTK('Property','Value',...) creates a new JTK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JTK_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JTK_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JTK

% Last Modified by GUIDE v2.5 18-Mar-2014 09:53:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JTK_OpeningFcn, ...
                   'gui_OutputFcn',  @JTK_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before JTK is made visible.
function JTK_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JTK (see VARARGIN)

% Choose default command line output for JTK
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JTK wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = JTK_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SelectFolder.
function SelectFolder_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% DataFolder=uigetdir('/Users/yishi/Documents/JTK/');
DataFolder=uigetdir('../');
set(handles.datafolder,'String',DataFolder);

% --- Executes on button press in AnalysisData.
function AnalysisData_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filefolder=get(handles.datafolder,'String');
filepath=strcat(filefolder,'/*.csv');
files=dir(filepath);
sz = size(files);

average_traces=cell(sz(1)+1);  %% store average traces, used for heatmap and dendrogram

for i=1:sz(1)+1
    average_traces{i}=zeros(1,145);
end

%%%if 'Analysis' folder doesn't exist, create it%%
if (exist('../Analysis','dir')~=7)
    mkdir('../Analysis');
end

csv_fid=fopen('../Analysis/averagetraces.csv','wt');  %%%csv file for heatmap
csv_fid2=fopen('../Analysis/averagetraces_cluster.csv','wt'); %%%csv file for dendrogram

average_traces{1,1}='filename';
for k=2:145
    average_traces{1,k}=(k-2)/6;
end

temperature_threshold=28; %%temperature threshold, temperature value over this threshold is valid

result_fid=fopen('../Analysis/finalresult.csv','wt');  %%store the final analysis result. 
fprintf(result_fid,'%s,%s,%s,%s\n','filename','JTK_PERIOD','JTK_LAG','JTK_AMP');

for mainloop=1:sz(1)
    handle = figure(mainloop);
%     set(gcf,'Units', 'inches', 'Position', [1 1 11 7]);   %%for mac
    set(gcf,'Units', 'inches', 'Position', [1 1 8 5]);  %%for windows
    hold on;
    
    filename=files(mainloop).name;
    str_writeback=strcat('../Analysis/','data_',filename,'.txt');
    writeback=fopen(str_writeback,'w+');
    fprintf(writeback,'%s\t','filename');
    
    
   str = strcat(filefolder,'/',filename);
   fid=fopen(str);
   title=textscan(fid,'%s',41,'delimiter',',');
   data=textscan(fid,'%s %s %f','delimiter',',');
   fclose(fid);
   
   time_label=data{1};  %%store all the time labels in one .csv file
   temperature=data{3}(:);  %%store all the temperature values in one .csv file
   [file_length,width]=size(time_label);
   
   %%compute moving average%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   mvaverage=zeros(1,file_length);
   windowsize=6; %%the gap is 1 hour
   mvaverage(1:windowsize-1)=temperature(1:windowsize-1);
   
   for w=windowsize:file_length
      value_in_window=temperature(w-windowsize+1:w);
      mvaverage(w)=mean(value_in_window);
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %%get the date of the first record temperature value%%%%%%%% 
   date_hour0=regexp(time_label(1),'\s+','split');
   date0=date_hour0{1}(1);
   time0=date_hour0{1}(2);
          
   date_day0=regexp(date0,'/','split');
   month0=date_day0{1}(1);
   day0=date_day0{1}(2);
   year0=date_day0{1}(3);
   firstmonth=str2double(month0{1});
   firstday=str2double(day0{1});
   firstyear=str2double(year0{1});
   firstyear=firstyear+2000;

   switch firstmonth
       case {2}
           
           if(isleapyear(firstyear))
               num_in_month=29;
           else
               num_in_month=28;
           end
           
           
           
       case {1,3,5,7,8,10,12}
           num_in_month=31;
           
       case {4,6,9,11}   
           num_in_month=30;
   end
   
   %%get the date of the last record temperature value%%%%%%%%
   %%then can compute  the total days the device recorded%%%%%
   date_hourN=regexp(time_label(file_length),'\s+','split');
   dateN=date_hourN{1}(1);
   timeN=date_hourN{1}(2);
          
   date_dayN=regexp(dateN,'/','split');
   monthN=date_dayN{1}(1);
   dayN=date_dayN{1}(2);
   lastmonth=str2double(monthN{1});
   lastday=str2double(dayN{1});
   
   if(firstmonth==lastmonth)
       day_index=1+lastday-firstday;  %%day_index stores the total days the device recorded
   else
       day_index=1+num_in_month-firstday+lastday;
   end
   
   %%%%get the first record that temperature value is over the threshold%%%
   for k=1:file_length
      if(temperature(k)>=temperature_threshold)
          date_hour1=regexp(time_label(k),'\s+','split');
          date1=date_hour1{1}(1);
          time1=date_hour1{1}(2);
          
          date_day1=regexp(date1,'/','split');
          day1=date_day1{1}(2);
 
          hour_minute1=regexp(time1,':','split');
          temphour1=hour_minute1{1}(1);
          tempminute1=hour_minute1{1}(2);
          hour1=str2double(temphour1{1});
          minute1=str2double(tempminute1{1});
          
          time_point1=hour1+minute1/60;
          breakpoint=k;
          break;
      end
   end
   
   %%%%%write all the valid value into 'data_filename.txt' file for analysis
   header=zeros();
   content=zeros();
   header(1)=time_point1;
   content(1)=temperature(breakpoint);
   hi=1;
   for k=(breakpoint+1):file_length
       if(temperature(k)>=temperature_threshold)
           hi=hi+1;
           header(hi)=header(hi-1)+1/6;
           content(hi)=temperature(k);
           
       else
           continue;
       end
   end
   
   for m=1:(hi-1)
       fprintf(writeback,'%f\t',header(m));
   end
   fprintf(writeback,'%f\r\n',header(hi));
   
   fprintf(writeback,'%s\t',filename);
   for n=1:hi
       fprintf(writeback,'%f\t',content(n));
   end
   
   fclose(writeback);
   
   columnnum=size(header);
   data_jtk=zeros(2,columnnum(2));
   data_jtk(1,:)=header;
   data_jtk(2,:)=content;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      for m=1:144
         if((time_point1-(m/6))<0)
             order=m;
             break;
         end
      end
      start_point=time_point1-(order-1)/6.0;
      
      time_index=zeros(1,day_index*144);
      for n=0:143
          time_index(n+1)=start_point + n/6;
      end
      
      reps=zeros(1,file_length);  %%the same as reps in R script, 
      
      temperature_each_day=zeros();
      
      for j=1:file_length
          if(temperature(j)>=temperature_threshold)
              
            reps(j)=1;
            date_hour=regexp(time_label(j),'\s+','split');
            date=date_hour{1}(1);
            time=date_hour{1}(2);
          
            hour_minute=regexp(time,':','split');
            temphour=hour_minute{1}(1);
            tempminute=hour_minute{1}(2);
            hour=str2double(temphour{1});
            minute=str2double(tempminute{1});
          
            time_point=hour+minute/60;
            point_index=int32(1+(time_point-start_point)*6);
 
            temperature_each_day(day_index,point_index)=temperature(j);
            
          else
              continue;
          end

      end
 
      %%compute the average temperature for one day 
      average_temperature=zeros(1,day_index*144);
      for jj=1:144
          temp_day_index=day_index;
          sum_temperature=0;
          for iii=1:day_index
              if(temperature_each_day(iii,jj)>0)
                  sum_temperature=sum_temperature+temperature_each_day(iii,jj);
              %else(temperature_each_day(iii,jj)==0)
              else
                  temp_day_index=temp_day_index-1;
              end
          end
           average_temperature(:,jj)=sum_temperature/temp_day_index;
      end
      
      %%repeat the average temperature for all the days%%%%%%%%%
      
      for ix=2:day_index
         average_temperature(1,(1+(ix-1)*144):144+(ix-1)*144)=average_temperature(1,1:144);
         time_index(1,(1+(ix-1)*144):144+(ix-1)*144)=time_index(1,1:144)+(ix-1)*24;
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      average_traces{(mainloop+1),1}=filename;
      for k=2:145
          average_traces{(mainloop+1),k}=average_temperature(1,k-1);
      end


      tempdate_hour=regexp(time_label(1),'\s+','split');
      tempdate=tempdate_hour{1}(1);
      temptime=tempdate_hour{1}(2);
      hour_min=regexp(temptime,':','split');
      temphour2=hour_min{1}(1);
      tempminute2=hour_min{1}(2);
      hour2=str2double(temphour2{1});
      minute2=str2double(tempminute2{1});
      first_point=hour2+minute2/60;
      pre_gap=int32(1+(first_point-start_point)*6);

      average_for_comparison=zeros(1,file_length);
      average_for_comparison(1:file_length)=average_temperature(pre_gap:(pre_gap+file_length-1));
      difference_avg_mvavg=abs(average_for_comparison-mvaverage);
  
      temperature2=temperature';
      
      time_index2=zeros(1,file_length);
      time_index2(1,1:file_length)=time_index(1,pre_gap:(pre_gap+file_length-1));
      
      
      %%%get the moving average value that is over the threshold%%
      zz=1;
      for iiii=1:file_length
          if(mvaverage(iiii)>=temperature_threshold)
              time_index3(zz)=time_index2(iiii);
              mvaverage2(zz)=mvaverage(iiii);
              zz=zz+1;
          else
              continue;
          end
      end
      
      time_index4=zeros(1,zz-1);
      mvaverage3=zeros(1,zz-1);
      time_index4(1:zz-1)=time_index3(1:zz-1);
      mvaverage3(1:zz-1)=mvaverage2(1:zz-1); 
      
%       avg=plot(time_index,average_temperature, 'Color', [0 0 1], 'LineWidth', 2); %plot average temperature
      avg=plot(time_index,average_temperature,'-b');
      raw=plot(time_index2,temperature2,'-r');
      mvavg=plot(time_index4,mvaverage3,'*g');
      
      legend([avg,raw,mvavg],'average','raw data','moving average',4);
      
    % Save the figure as a .tif file at 300 DPI.
    
      fig_name=strcat('../Analysis/',filename);
      print(handle, [fig_name, '.tif'],'-dtiff','-r300');
      hold off;  
    
    
    
    %%%%%%%%plot moving average after hpfilter%%%%%%%%%%%%%
    
%     tempindex=size(time_index4);
%     xdata=zeros(tempindex(2));
%     xdata_raw=time_index4;
%     ydata=zeros(tempindex(2));
%     ydata=mvaverage3;
%     xdata=xdata_raw;
%     [ydata2] = hpfilter(ydata,1600);
% 
%     handle2 = figure(mainloop+sz(1));
%     hold on;
%     hp=plot(xdata,ydata2,'-g');
%     legend(hp,'moving average after hp filter',4);
%     fig_name2=strcat('../Analysis/',filename,'_hpfilter');
%     print(handle2, [fig_name2, '.tif'],'-dtiff','-r300');
%     hold off;

    %%%%%%%%end for hpfilter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%use JTK_CYCLE to compute lag, amplitude and period%%%%%%%    
    
    data_length=sum(reps);
    fid_data=fopen(str_writeback,'rt');
    firstrow_name=fscanf(fid_data,'%s',1);
    firstrow_data=fscanf(fid_data,'%f',data_length);
    secondrow_name=fscanf(fid_data,'%s',1);
    secondrow_data=fscanf(fid_data,'%f',data_length);
    fclose(fid_data);
    
    
    fprintf(result_fid,'%s,',filename);
    
    JTK_AMPFACTOR=sqrt(2);
    JTK_PIHAT= 3.1416;
    timepoints=file_length;
    
    %%%%jtkdist, jtk_init, jtkx functions are just the same as in R script%
    [JTK_GRP_SIZE,JTK_NUM_GRPS,JTK_NUM_VALS,JTK_MAX,JTK_DIMS,JTK_SDV,JTK_EXV,JTK_EXACT,JTK_CP]=jtkdist(timepoints,reps);
    periods=144:144;
    interval=0.16667;
    [JTK_CGOOSV,JTK_SIGNCOS,JTK_PERIODS,JTK_INTERVAL]=jtk_init(periods,interval,JTK_GRP_SIZE,JTK_NUM_GRPS,JTK_NUM_VALS,JTK_DIMS,JTK_PIHAT);
    
    conf=0.9;
    z=secondrow_data;
    [JTK_PERIOD,JTK_LAG,JTK_AMP]=jtkx(z,conf,JTK_DIMS,JTK_PERIODS,JTK_CGOOSV,JTK_SIGNCOS,JTK_MAX,JTK_EXACT,JTK_EXV,JTK_SDV,JTK_INTERVAL,JTK_AMPFACTOR,JTK_CP);
  
    fprintf(result_fid,'%f,%f,%f\n',JTK_PERIOD,JTK_LAG,JTK_AMP);
    
    
    %%%%%%%%%%%%end for compute amplitude, lag, and period%%%%%%%%%%%%%%%%%%
    
    
    
end  %%%end for the main loop,   
   
fclose(result_fid);

%%%write the average traces into related files used for plot heatmap and
%%%dendrogram

for i=1:sz(1)+1
    for j=1:144
        avt=num2str(average_traces{i,j});
        fprintf(csv_fid,'%s,',avt);
    end
    avt=num2str(average_traces{i,145});
    fprintf(csv_fid,'%s\n',avt);
    
end


for i=2:sz(1)+1
    for j=1:144
        avt2=num2str(average_traces{i,j});
        fprintf(csv_fid2,'%s,',avt2);
    end
    avt2=num2str(average_traces{i,145});
    fprintf(csv_fid2,'%s\n',avt2);
    
end


fclose(csv_fid);
fclose(csv_fid2);

display('Analysis finished, you can run heatmap and dendrogram.');  %%end for 'Analysis' button





function datafolder_Callback(hObject, eventdata, handles)
% hObject    handle to datafolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datafolder as text
%        str2double(get(hObject,'String')) returns contents of datafolder as a double


% --- Executes during object creation, after setting all properties.
function datafolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datafolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in heatmap.
function heatmap_Callback(hObject, eventdata, handles)
% hObject    handle to heatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Creates a heatmap for a given set of data, ranging from blue (low) to
% yellow (high), with black in between.  
%
% *** Note-- you'll probably have to change the location of the box around
% the legend if you have anything other than 13 columns.  I can make it
% more automatic, or you can comment out the rectangle (line 135) and 
% draw it in yourself.
%
% Inputs: 
% - data: name of a .csv file that contains the data, with row/column
% names.  It'll object if you don't have headers.
% - Min: the minimum value to plot.  If the score is less than the minimum
% value, it'll appear blue.
% - Max: the maximum value to plot.  If the score is greater than the 
% maximum, it'll appear yellow.
%
% A score halfway between the min and max will appear black.
%


data_source='../Analysis/averagetraces.csv';


Min=0.5;   % two values that you need to set
Max=1.5;   % 


% Import data.
allData_heatmap = importdata(data_source,',');

 data_heatmap = allData_heatmap.data; % Numeric data
 rows_heatmap = allData_heatmap.textdata(2:size(allData_heatmap.textdata,1),1); % Row names
 cols_heatmap = allData_heatmap.textdata(1,2:size(allData_heatmap.textdata,2)); % Column names
 
 if ~ length(cols_heatmap) % If the first row of the data is numerical rather than text, it'll be imported in as data-- so remove it.
     data_heatmap(1,:) = [];
 end
 
height = size(data_heatmap,1); % # rows
width = size(data_heatmap, 2); % # cols

median_val=0;
for i=1:height
    median_val=median(data_heatmap(i,:));
    for j=1:width
        data_heatmap(i,j)=data_heatmap(i,j)/median_val;
    end
end




boxSize = 15;

%figure('Position', [1     1   464   597]);

%set(gca,'Position', [0.01 0.01 .99 .99]);
figure();

% set(gcf,'Units', 'inches', 'Position', [1 1 11 7]); %%for mac
set(gcf,'Units', 'inches', 'Position', [1 1 8 5]);  %%for windows
hold on; % Allows two graphs to be printed on top of each other.


for i=1:height
    for j = 1:width
        score = (data_heatmap(i,j)-Min)/(Max-Min);  % Normalize so the scores range from 0 to 1.
        
        if (score > 1) % If it's greater than the Max value, set to the max.
            score = 1;
        end

        if (score < 0) % If it's less than the Min value, set to the min.
            score = 0;
        end

        % Determine the color to plot, in RGB.  
        % If the score is zero (i.e. less than the minimum), it'll appear blue (RGB = [0 0 1]).
        % If the score is 0.5 (halfway b/w min and max), it'll appear black (RGB = [0 0 0]).
        % If the score is one (i.e. greater than the maximum), it'll appear yellow (RGB = [1 1 0])
        
        % Between those three cases, the color will be some shade of
        % blue-black (RGB = [0 0 b], where 0 <= b <= 1, b = 1 - 2*score);
        % or yellow-black (RGB = [y y 0], where 0 <= y <= 1, y = 2*score - 1);
        
        if (score <= 0.5)
            R = 0;
            G = 0;
            B = -2*score+1;
        else
            R = 2*score-1;
            G = 2*score - 1;
            B = 0;
        end
        
            color = [R G B];
        
        % Plot each point as a square.
       % plot(j*2,-i*2, 's', 'MarkerSize', 15, 'MarkerEdgeColor', color, 'MarkerFaceColor', color)
       rectangle('Position', [(j-1)*boxSize, -(i-1)*boxSize, boxSize, boxSize], 'FaceColor', color, 'EdgeColor', 'none');
    end
end

colNames{1} = '';
rowNames{1} = '';

for k = 1:width
%    colNames{k+1} = cols{k};
end

counter = 2;
for l = height:-1:1
 %   rowNames{counter} = rows{l};
  %  counter = counter + 1;
end


% Plot properties so it looks pretty.
% Resize axes.
axis equal; % Forces each data point to be square.  Comment this out if you want to be able to make the points rectangular as you resize the plot.

axis([0 (width)*boxSize -boxSize*(height-1) boxSize])

% Resize tick marks (one per datapoint).
%set(gca, 'XTick', [0:2:2*(width)]);
%set(gca, 'YTick', [-2*(height+1):2:0]);

% Font size of axis labels.
set(gca,'FontSize',8)

% Label x & y axes with the row/column names.
%set(gca,'XTickLabel',colNames)
 %xticklabel_rotate([],45,[],'Fontsize',8)
%set(gca,'YTickLabel',rowNames)


% Plot scale bar.

% plot(2*(width+1)+2, -2, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1]); %Score = 0.000.
% text(2*(width+1)+3, -2, ['\leq ', num2str(Min)])
% plot(2*(width+1)+2, -3, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0.875], 'MarkerFaceColor', [0 0 0.875]);
% plot(2*(width+1)+2, -4, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0.75], 'MarkerFaceColor', [0 0 0.75]); %Score = 0.125.
% plot(2*(width+1)+2, -5, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0.625], 'MarkerFaceColor', [0 0 0.625]);
% plot(2*(width+1)+2, -6, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0.5], 'MarkerFaceColor', [0 0 0.5]); %Score = 0.250.
% text(2*(width+1)+3, -6, num2str(0.25*(Max-Min)+Min))
% plot(2*(width+1)+2, -7, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0.375], 'MarkerFaceColor', [0 0 0.375]); 
% plot(2*(width+1)+2, -8, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0.25], 'MarkerFaceColor', [0 0 0.25]); %Score = 0.375.
% plot(2*(width+1)+2, -9, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0.125], 'MarkerFaceColor', [0 0 0.125]);
% plot(2*(width+1)+2, -10, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]); %Score = 0.500.
% text(2*(width+1)+3, -10, num2str(0.5*(Max-Min)+Min))
% plot(2*(width+1)+2, -11, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0.125 .125 0], 'MarkerFaceColor', [0.125 .125 0]); 
% plot(2*(width+1)+2, -12, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0.25 .25 0], 'MarkerFaceColor', [0.25 .25 0]); %Score = 0.625.
% plot(2*(width+1)+2, -13, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0.375 .375 0], 'MarkerFaceColor', [0.375 .375 0]); 
% plot(2*(width+1)+2, -14, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0.5 .5 0], 'MarkerFaceColor', [.5 .5 0]); %Score = 0.750.
% text(2*(width+1)+3, -14, num2str(0.75*(Max-Min)+Min))
% plot(2*(width+1)+2, -15, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0.625 .625 0], 'MarkerFaceColor', [0.625 .625 0]); 
% plot(2*(width+1)+2, -16, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [.75 .75 0], 'MarkerFaceColor', [.75 .75 0]); %Score = 0.875.
% plot(2*(width+1)+2, -17, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [0.875 .875 0], 'MarkerFaceColor', [0.875 .875 0]); 
% plot(2*(width+1)+2, -18, 's', 'MarkerSize', 25, 'MarkerEdgeColor', [1 1 0], 'MarkerFaceColor', [1 1 0]); %Score = 1.000.
% text(2*(width+1)+3, -18, ['\geq ', num2str(Max)])

%rectangle('Position', [2*(width+1)+1,-2*(height)+7, 4, height+5]) % Draws a box around the scale bar.

set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gca, 'XTickLabel',[]);
set(gca, 'YTickLabel', []);
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');

hold off;





% --- Executes on button press in dendrogram.
function dendrogram_Callback(hObject, eventdata, handles)
% hObject    handle to dendrogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Have the user select a .csv file.
%[file, path, ~] = uigetfile('*.csv', 'Select the .csv file', 'MultiSelect', 'off');

%Import data.
%allData = importdata([path, file],',');


data_cluster_source='../Analysis/averagetraces_cluster.csv';
%data='/Users/yishi/Documents/test/averaged_traces_cluster.csv';

% Import data.
allData_cluster = importdata(data_cluster_source,',');

% Fish out the data.
data_cluster = allData_cluster.data;

Z_cluster = linkage(data_cluster, 'ward');
figure();
hd=dendrogram(Z_cluster, 0, 'Labels', allData_cluster.rowheaders, 'orientation', 'left','ColorThreshold','default');
set(hd,'LineWidth',2);
