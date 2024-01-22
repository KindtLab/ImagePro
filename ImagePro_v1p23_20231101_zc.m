function varargout = ImagePro_v1p23_20231101_zc(varargin)
% IMAGEPRO_V1P23_20231101_ZC M-file for ImagePro_v1p23_20231101_zc.fig
%      IMAGEPRO_V1P23_20231101_ZC, by itself, creates a new IMAGEPRO_V1P23_20231101_ZC or raises the existing
%      singleton*.
%
%      H = IMAGEPRO_V1P23_20231101_ZC returns the handle to a new IMAGEPRO_V1P23_20231101_ZC or the handle to
%      the existing singleton*.
%
%      IMAGEPRO_V1P23_20231101_ZC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEPRO_V1P23_20231101_ZC.M with the given input arguments.
%
%      IMAGEPRO_V1P23_20231101_ZC('Property','Value',...) creates a new IMAGEPRO_V1P23_20231101_ZC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImagePro_v1p23_20231101_zc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImagePro_v1p23_20231101_zc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help ImagePro_v1p23_20231101_zc
% Last Modified by GUIDE v2.5 01-Nov-2023 12:53:23
% Begin initialization code - DO NOT EDIT

% Update Log:
% leiz 20231101
% fixed Movie button for saving current images in the listbox into .mp4
% fixed colorbar button



gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImagePro_v1p23_20231101_zc_OpeningFcn, ...
                   'gui_OutputFcn',  @ImagePro_v1p23_20231101_zc_OutputFcn, ...
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
% --- Executes just before ImagePro_v1p23_20231101_zc is made visible.
function ImagePro_v1p23_20231101_zc_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to ImagePro_v1p23_20231101_zc (see VARARGIN)
% Choose default command line output for ImagePro_v1p23_20231101_zc
handles.output = hObject;
%
colorTmp=[1,0,0;
          0,1,0;
          0,0,1;
          0.502,0,1;
          0,1,1;
          1,0,1;
          0.3,0.3,0.7;
          0.5,1/2,1/2];


handles.appBG = [];
handles.colorTmp=colorTmp;
handles.filePathMat={};
handles.fileNameMat={};
handles.fileMat=[];
handles.pMat=[];
handles.appHome = fileparts(which(mfilename));
handles.bgColor=get(gcf,'color');

if exist([handles.appHome,filesep,'zfhc.png'])
    axis(handles.axes_rawImg);
    imshow(imread([handles.appHome,filesep,'zfhc.png']));
else
    text(0.5, 0.6, {'\bfKindt Lab';'\fontsize{36}For Fish Fun'},...
        'interpreter','tex','Units','normalized',...
        'horizontalAlignment','center','Fontsize',52,'Color',[0.75,0.75,0.75]);
end

handles.fileType = 'png';

axis off;
hold off;
guidata(hObject, handles);
% UIWAIT makes ImagePro_v1p23_20231101_zc wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = ImagePro_v1p23_20231101_zc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;
% --- Executes on button press in pb_imgLoad.
function pb_imgLoad_Callback(hObject, eventdata, handles)
% hObject    handle to pb_imgLoad (see GCBO)
 cdDefault=cd;
  set(handles.pb_imgLoad,'backgroundcolor',handles.bgColor);
 % t=get(gcf,'children');
 %  set(findobj(t,'style','text'),'backgroundcolor',c1);
 %   set(findobj(t,'style','checkbox'),'backgroundcolor',c1);
 % set(findobj(t,'style','pushbutton'),'backgroundcolor',c1);
 % set(findobj(t,'type','uipanel'),'backgroundcolor',c1)
method='imgLoad';
handles=imgLoadInitial(hObject,handles,method);
if strcmp(handles.fileType,'mat')~=1
    cd(handles.filePath)
    imgInfo=imfinfo(handles.fileName(1).name);
    uint8Flag=imgInfo.BitDepth;
    if uint8Flag==24
        uint8Flag=8;
    elseif uint8Flag==16*3
        uint8Flag=16;
    elseif uint8Flag==1
        uint8Flag=8;
    end
    handles.imgDepth=['uint',num2str(uint8Flag)];
else
    handles.imgDepth='uint8';
end
handles=pp_spfilterMethod_Callback(handles.pp_spfilterMethod, eventdata, handles);
% handles=pp_2filterMethod_outside(hObject,handles);
% handles=filterInfoRead2_spatial(handles);
% cd(cdDefault);
%set(allchild(gcf),'unit','normalized')
set(handles.pb_imgLoad,'BackgroundColor',[0.3,0.75,0.93]);
guidata(hObject,handles);

function handles=imgLoadInitial(hObject,handles,method)
pathWhich=which(mfilename);
dotNO=strfind(pathWhich,filesep);
dotNO=dotNO(end);
pathWhich=pathWhich(1:dotNO-1);
path(pathWhich,path);
cdDef=cd;

%% reset .mat file list
handles.filePathMat={};
handles.fileNameMat={};
handles.fileMat=[];

%%
if strcmp(method,'imgLoad')
    if isfield(handles,'filePath_output')
        cd(handles.filePath_output)
    end
    if get(handles.cb_cdCurrent,'value')
        cd(cdDef)
    end
    [file,filePath]=uigetfile('*.*');
    if filePath==0
        filePath=handles.filePath_output;
        file=handles.file;
    end
    cd(filePath);
    indexNO=strfind(file,'.');
    fileType=file(indexNO(end)+1:end);
    fileName=dir(['*.',fileType]);
    for iss=length(fileName):-1:1
        if strcmp(fileName(iss).name(1:2),'._')
            fileName(iss)=[];
        end
    end
    p=length(fileName);
elseif strcmp(method,'userInput')
    filePath=get(handles.et_fileName_path,'string');
    file1=get(handles.et_fileName_file,'string');
    file2=get(handles.et_fileName_start,'string');
    file3=get(handles.et_fileName_format,'string');
    % check whether file3 contains dot
    indexNO=strfind(file3,'.');
    if isempty(indexNO)
        file3=['.',file3];
    end
    fileType=file3(2:end);
    file=[file1,file2,file3];
    p=eval(get(handles.et_fileName_to,'string'))-eval(file2)+1;
    zeroPadding=length(file2);
    fileName=struct;
    startNO=eval(file2);
    step=eval(get(handles.et_step,'string'));
    ss=1;
    if get(handles.cb_NoZero,'value')
    for ii=1:step:p
        jj=startNO+ii-1;
        file2=num2str(jj);
        fileName(ss).name=[file1,file2,file3];
        ss=ss+1;
    end
    else
         for ii=1:step:p
            jj=startNO+ii-1;
            file2=num2str(jj,['%0',num2str(zeroPadding),'d']);
            fileName(ss).name=[file1,file2,file3];
            ss=ss+1;
        end
    end
    file=fileName(1).name;
end
p=length(fileName);
%%
%show images' name
handles.p=p;
handles.fileType=fileType;
handles.filePath_output = filePath;
handles.filePathOri = filePath;
handles.file=file;
handles.filePath=filePath;
handles.fileName=fileName;
handles.filePath_movie = [filePath,filesep,'Movies'];
switch fileType
    case ('mat')
        handles.filePathMat=filePath;
        handles.fileNameMat=fileName;
        handles.fileMat=fileName(1).name;
        handles.filePath_output = fileparts(fileparts(filePath));
end

listName=cell(p,1);
for ii=1:p
    listName{ii}=fileName(ii).name;
    
end
set(handles.lb_curFileNames,'value',1);
set(handles.lb_curFileNames,'string',listName);
%  filePathWrap=textwrap(handles.st_curFilePath,{filePath});
filePathWrap=filePath;
set(handles.st_curFilePath,'string',filePathWrap);
set(handles.t_curFile,'string',file)
handles=imageShow(handles,'replace');
set(handles.et_IOSresultPath,'string',handles.filePath_output);
contents = cellstr(get(handles.pm_groupNO,'String'));
if isfield(handles,'pointPst')==0
    handles.pointPst=cell(str2num(contents{get(handles.pm_groupNO,'Value')}),1);
end
set(handles.dt_IntestedImg,'string',['[1:',num2str(p),']'])
%% for showing point coordinateions
set(allchild(gca),'ButtonDownFcn', ...
    'IOS_Software(''axes_rawImg_ButtonDownFcn'',gco,[],guidata(gcbo))')
set(allchild(gca),'hitTest','on')
%% for spatial binning
cd(handles.filePath)
img=differentTypeRead(handles.file,handles.fileType);
[mm,nn,nnTrash]=size(img);
% set(handles.et_xBinningEnd,'string',num2str(nn))
% set(handles.et_yBinningEnd,'string',num2str(mm))
function img=differentTypeRead(file,fileType)
if strcmp(fileType,'txt')
    img=load(file);
elseif strcmp(fileType,'mat')
    img2=load(file);
    img2Name=fieldnames(img2);
    img=img2.(img2Name{1});
else
    img=single(imread(file));
%     img=double(imread(file));
%     if size(img,3)~=1
%         img=img(:,:,1);
%     end
end
function img=differentTypeReadFilter(img,filterMethod2)
% if strcmp(fileType,'txt')
%     img=load(file);
% elseif strcmp(fileType,'mat')
%     img2=load(file);
%     img2Name=fieldnames(img2);
%     img=img2.(img2Name{1});
% else
%     img=single(imread(file));
% %     img=double(imread(file));
% %     if size(img,3)~=1
% %         img=img(:,:,1);
% %     end
% end
sizNo=filterMethod2.sizNo;
sigmaNo=filterMethod2.sigmaNo;
nameNo=filterMethod2.name;
%% filterImg
if nameNo==0
    %imgAvg=imgAvg;
elseif nameNo==1
    w=fspecial('gaussian',[sizNo,sizNo],sigmaNo);
    img=imfilter(img(:,:,1),w,'replicate');
elseif nameNo==2
    img=medfilt2(img(:,:,1),[sizNo,sizNo]);
elseif nameNo==3
        img=wiener2(img(:,:,1),[sizNo,sizNo]);
elseif nameNo==4
    imgAvg2=maxFilter2(img(:,:,1),sizNo);
    img=single(imgAvg2);
%     imgAvg2=edge(img(:,:,1),'canny',[.05,.4],sigmaNo);
%     img=single(imgAvg2);
elseif nameNo==5
   order=filterMethod2.order;
   fcut=filterMethod2.fcut;
   img=butter_fcut(img,fcut,1,order,'abs');
end
function [zoFilter,varargout]=butter_fcut(z0,fcut,method,butterN,rc)
% method=1 butterworth;
[m,n]=size(z0);
mm=1:m;
nn=1:n;
m0=m/2;
if mod(m0,2)==0
    m0=m0+0.5;
else
    m0=floor(m0)+1;
end

n0=n/2;
if mod(n0,2)==0
    n0=n0+0.5;
else
    n0=floor(n0)+1;
end
[xx,yy]=meshgrid(nn,mm);
yy=yy-m0;
xx=xx-n0;
dist=sqrt(yy.^2+xx.^2);
if method==1

    map=1 ./ (1.0 + (dist ./ (fcut*m)).^(2*butterN));
    zoFilter=ifft2(ifftshift(map.*fftshift(fft2(z0))));
elseif method==2
    map=dist<fcut*m;
    zoFilter=ifft2(ifftshift(map.*fftshift(fft2(z0))));    
end
if strcmp(rc,'abs')
    zoFilter=abs(zoFilter);
end
if nargout>=2
    varargout{1}=map;
end
function  img2=maxFilter2(img1,sizNo)
[mm,nn]=size(img1);
sizNo=floor(sizNo/2);
edgeY=sizNo;
edgeX=sizNo;
%% make the images bigger; replicate padding
% img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
% img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
img1_big(1:edgeY,:)=ones(edgeY,1)*img1_big(1+edgeY,:);
img1_big(edgeY+mm+1:end,:)=ones(edgeY,1)*img1_big(mm,:);
img1_big(:,1:edgeX)=img1_big(:,edgeX+1)*ones(1,edgeX);
img1_big(:,nn+1+edgeX:end)=img1_big(:,nn)*ones(1,edgeX);
sizeY1=edgeY;
sizeX1=edgeX;
%% initialize matrix
img2=zeros(mm,nn);
for iiTmp=1:mm
%          waitbar(iiTmp/mm,h_wait,[num2stimageShowr(100*iiTmp/mm,'%04.1f'),'%completed']);
    ii=edgeY+iiTmp;
    for jjTmp=1:nn
        jj=edgeX+jjTmp;
        ROI=img1_big(ii-sizeY1:ii+sizeY1,jj-sizeX1:jj+sizeX1); %#ok<*PFBNS>
        img2(iiTmp,jjTmp)=max(ROI(:));
%         disp(['y',num2str(iiTmp,'%03d'),'jj',num2str(jjTmp,'%03d')])
    end
end
function handles=imageShow(handles,varargin)
contents = cellstr(get(handles.pm_colormap_raw,'String'));
colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
if get(handles.cb_uint16,'value')
colorSelectedTmp{1}=[colorSelectedTmp{1}(1:end-5),'(2^16)'];
end
colorSelected=eval(colorSelectedTmp{1});
% cd(handles.filePath)
handles.newIntensity=inline(get(handles.et_inline,'string'));
fileType=handles.fileType;
file=handles.file;
img=differentTypeRead(fullfile(handles.filePath,file),fileType);
handles.imgRaw=img;
 newImg=handles.newIntensity(img);


%% further processing
if get(handles.cb_eval,'value')
    x=newImg;
    eval(get(handles.et_eval,'string')) ;
    newImg=x;
end
handles.img=newImg;
d3=size(newImg,3);
if nargin==2 && strcmp(varargin{1},'replace')
    axes(handles.axes_rawImg);
    % set(gcf,'currentAxes',handles.axes_rawImg)
    hold off;
    if d3==3
        if get(handles.cb_uint16,'value')
            imshow(uint16(newImg));
            newImgColor=uint16(newImg);
        else
            imshow(uint8(newImg))
            newImgColor=uint8(newImg);
        end
        
    elseif d3==1

       %imshow(uint8(newImg),colorSelected)
        if get(handles.cb_uint16,'value')
            newImgColor=ind2rgb(uint16(newImg),colorSelected);
            imshow(newImgColor)
        else
            newImgColor=ind2rgb(uint8(newImg),colorSelected);
            imshow(newImgColor)
        end       
    end
else
       h_img=findobj(handles.axes_rawImg,'type','image');
    %set(h_img,'CData',uint16(img))
    if d3==3
            if get(handles.cb_uint16,'value')
                set(h_img,'CData',uint16(newImg))
                newImgColor=uint16(newImg);
            else
                set(h_img,'CData',uint8(newImg))
                newImgColor=uint8(newImg);
            end        
         
    elseif d3==1
        newImgColor=ind2rgb(uint8(newImg),colorSelected);
       set(h_img,'CData',newImgColor)
            if get(handles.cb_uint16,'value')
                newImgColor=ind2rgb(uint16(newImg),colorSelected);
               set(h_img,'CData',newImgColor)
            else
            newImgColor=ind2rgb(uint8(newImg),colorSelected);
           set(h_img,'CData',newImgColor)
            end          
      %  set(h_img,'CData',uint8(newImg))
    end
end

% axes(handles.axes_rawImg);
% colormap(colorSelected);
% colorbar;

img=newImg;
img=double(img);
[m,n,nn_trash]=size(img);
imgQuality=cell(5,1);
imgQuality{1}=[num2str(m),'x',num2str(n)];
imgQuality{2}=num2str(mean(img(:)));
imgQuality{3}=num2str(std(img(:)));
imgQuality{4}=num2str(max(img(:)));
imgQuality{5}=num2str(min(img(:)));
set(handles.imgQuality,'data',imgQuality)
handles.newImgColor=newImgColor;
% --- Executes on selection change in lb_curFileNames.
function handles=lb_curFileNames_Callback(hObject, eventdata, handles)
% hObject    handle to lb_curFileNames (see GCBO)
contents = get(hObject,'String');
if length(contents)<get(hObject,'Value')
    set(hObject,'value',1)
end
file=contents{get(hObject,'Value')};
handles.file=file;
handles=imageShow(handles);
set(handles.t_curFile,'string',file);
set(handles.et_file,'string',num2str(get(hObject,'Value')));
if isfield(handles,'pointer')
    if ~isempty(handles.pointer)
        handles=IOS_time_plotCurve(hObject,handles);
    end
end
guidata(hObject,handles);
function lb_curFileNames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_curFileNames (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: listbox controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_file_Callback(hObject, eventdata, handles)
% hObject    handle to et_file (see GCBO)
fileNO=str2double(get(hObject,'String'));
set(handles.lb_curFileNames,'value',fileNO);
lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles);
% et_file as text
%        str2double(get(hObject,'String')) returns contents of et_file as a double
function et_file_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to et_file (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_inline_Callback(hObject, eventdata, handles)
% hObject    handle to et_inline (see GCBO)
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles)
% et_inline as text
%        str2double(get(hObject,'String')) returns contents of et_inline as a double
function et_inline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_inline (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_IOSresultPath_Callback(hObject, eventdata, handles)
function et_IOSresultPath_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_pre1_Callback(hObject, eventdata, handles)
function et_pre1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_pre2_Callback(hObject, eventdata, handles)
function et_pre2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%
function handles=IOS_time_parameterRead(handles)
handles.IOSresultPath=get(handles.et_IOSresultPath,'string');
handles.pre1=round(eval(get(handles.et_pre1,'string')));
handles.pre2=round(eval(get(handles.et_pre2,'string')));
% handles.x1=eval(get(handles.et_x1,'string'));
% handles.x2=eval(get(handles.et_x2,'string'));
% handles.NO=round(eval(get(handles.et_NO,'string')));
% handles.method=get(handles.pm_method,'value');
function et_rawIOSfolderName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_rawIOS_mat.
function handles=pb_rawIOS_mat_Callback(hObject, eventdata, handles)
handles=IOS_time_parameterRead(handles);
pre1=handles.pre1;
pre2=handles.pre2;
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
fileName=handles.fileName;
cd(filePath);
img=differentTypeRead(file,fileType);
p=handles.p;
inverseNegative=get(handles.cb_inverseNegative,'value');
%% prestimulus mean
baselineImg=zeros(size(img));
for ii=pre1:pre2
    file=fileName(ii).name;
    img=differentTypeRead(file,fileType);
    baselineImg=baselineImg+img;
end
baselineImg=baselineImg/(pre2-pre1+1);
baselineImg(baselineImg==0)=1;
%% IOS
IOSresultPath=handles.filePath_output;
resultName=get(handles.et_rawIOSfolderName_mat,'string');
resultPath=fullfile(IOSresultPath,resultName);
if  exist(resultPath)~=7
    %     cd(IOSresultPath)
    mkdir(resultPath)
end
fileTypeLength=length(fileType);
hh=waitbar(0,'just wait','name','Saving dF/F0 in .mat');
binningN=eval(get(handles.et_binning,'string'));
pp=floor(p/binningN);
for ii=1:pp
    waitbar(ii/pp,hh, ...
        [num2str(100*ii/pp,'%03.1f'),'% completed']);
    cd(filePath)
    imgAvg=img*0;
    for jj=1:binningN
        file=fileName((ii-1)*binningN+jj).name;
        img=differentTypeRead(file,fileType);
        imgAvg=imgAvg+img;
    end
    imgAvg=imgAvg/binningN;
    if get(handles.cb_dF,'value')
        imgIOStmp=(imgAvg-baselineImg);
    else
        imgIOStmp=(imgAvg-baselineImg)./baselineImg;
    end
    
    %imgIOStmp=(imgAvg-baselineImg);
    if inverseNegative==1
        imgIOStmp=abs(imgIOStmp);
    end
    imgIOS=single(imgIOStmp);
    cd(resultPath)
    saveName=[file(1:end-fileTypeLength),'mat'];
    save(saveName,'imgIOS');
end
close(hh)
%% save them for other applications, i.e., creating grayscale & color images
handles.mat_resultPath=resultPath;
handles.mat_fileType='mat';
handles.mat_file=saveName;
guidata(hObject,handles)
function et_rawIOSfolderName_mat_Callback(hObject, eventdata, handles)
function et_rawIOSfolderName_mat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
function et_rawIOSfolderName_Callback(hObject, eventdata, handles)
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_g1.
function pb_g1_Callback(hObject, eventdata, handles)
gN=1;
handles=getPointG1(handles,gN);
guidata(hObject,handles);
% --- Executes on button press in pb_g2.
function pb_g2_Callback(hObject, eventdata, handles)
gN=2;
handles=getPointG1(handles,gN);
guidata(hObject,handles);
% --- Executes on button press in pb_g3.
function pb_g3_Callback(hObject, eventdata, handles)
gN=3;
handles=getPointG1(handles,gN);
guidata(hObject,handles);
% --- Executes on button press in pb_g4.
function pb_g4_Callback(hObject, eventdata, handles)
gN=4;
handles=getPointG1(handles,gN);
guidata(hObject,handles);
% --- Executes on button press in pb_g5.
function pb_g5_Callback(hObject, eventdata, handles)
gN=5;
handles=getPointG1(handles,gN);
guidata(hObject,handles);
function pushbutton10_Callback(hObject, eventdata, handles)
gN=handles.pointCurrent_group;
handles.pointPst{gN(end)}(end)=[];
delete(handles.pointCurrent_h(end));
handles.pointCurrent_h(end)=[];
handles.pointCurrent_group(end)=[];
handles.pointPstH{gN(end)}(end)=[];
handles=GroupROITag(handles);
guidata(hObject,handles)
% --- Executes on button press in pb_g6.
function pb_g6_Callback(hObject, eventdata, handles)
gN=6;
handles=getPointG1(handles,gN);
guidata(hObject,handles);
function pm_groupNO_Callback(hObject, eventdata, handles)
contents = cellstr(get(handles.pm_groupNO,'String'));
if get(handles.pm_groupNO,'value')<7
    groupNO=eval(contents{get(handles.pm_groupNO,'Value')});
    handles.pointPst=cell(groupNO,1);
    handles.pointPstH=cell(groupNO,1);
    if groupNO==1
        set(handles.pb_g1,'visible','on')
        set(handles.pb_g2,'visible','off')
        set(handles.pb_g3,'visible','off')
        set(handles.pb_g4,'visible','off')
        set(handles.pb_g5,'visible','off')
        set(handles.pb_g6,'visible','off')
        set(handles.pb_g7,'visible','off')
        set(handles.pb_g7_et,'visible','off')
    elseif groupNO==2
         set(handles.pb_g1,'visible','on')
        set(handles.pb_g2,'visible','on')
        set(handles.pb_g3,'visible','off')
        set(handles.pb_g4,'visible','off')
        set(handles.pb_g5,'visible','off')
        set(handles.pb_g6,'visible','off')
         set(handles.pb_g7,'visible','off')
        set(handles.pb_g7_et,'visible','off')
    elseif groupNO==3
         set(handles.pb_g1,'visible','on')
        set(handles.pb_g2,'visible','on')
        set(handles.pb_g3,'visible','on')
        set(handles.pb_g4,'visible','off')
        set(handles.pb_g5,'visible','off')
        set(handles.pb_g6,'visible','off')
         set(handles.pb_g7,'visible','off')
        set(handles.pb_g7_et,'visible','off')
    elseif groupNO==4
         set(handles.pb_g1,'visible','on')
        set(handles.pb_g2,'visible','on')
        set(handles.pb_g3,'visible','on')
        set(handles.pb_g4,'visible','on')
        set(handles.pb_g5,'visible','off')
        set(handles.pb_g6,'visible','off')
        set(handles.pb_g7,'visible','off')
        set(handles.pb_g7_et,'visible','off')
    elseif groupNO==5
         set(handles.pb_g1,'visible','on')
        set(handles.pb_g2,'visible','on')
        set(handles.pb_g3,'visible','on')
        set(handles.pb_g4,'visible','on')
        set(handles.pb_g5,'visible','on')
        set(handles.pb_g6,'visible','off')
        set(handles.pb_g7,'visible','off')
        set(handles.pb_g7_et,'visible','off')
    elseif groupNO==6
         set(handles.pb_g1,'visible','on')
        set(handles.pb_g2,'visible','on')
        set(handles.pb_g3,'visible','on')
        set(handles.pb_g4,'visible','on')
        set(handles.pb_g5,'visible','on')
        set(handles.pb_g6,'visible','on')
        set(handles.pb_g7,'visible','off')
        set(handles.pb_g7_et,'visible','off')
    end
elseif get(handles.pm_groupNO,'value')==7
        set(handles.pb_g1,'visible','off')
        set(handles.pb_g2,'visible','off')
        set(handles.pb_g3,'visible','off')
        set(handles.pb_g4,'visible','off')
        set(handles.pb_g5,'visible','off')
        set(handles.pb_g6,'visible','off')
        set(handles.pb_g7,'visible','on')
        set(handles.pb_g7_et,'visible','on')
        handles.pointPst=cell(7,1);
        handles.pointPstH=cell(7,1);
end
guidata(hObject,handles)
function pm_groupNO_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function handles=getPointG1(handles,gN)
colorTmp=handles.colorTmp;
pointPst=handles.pointPst;
 axes(handles.axes_rawImg);
% set(gcf,'currentAxes',handles.axes_rawImg)
% zoom on
% h=zoom;
% set(h,'ActionPostCallback','uiresume(gcbf)');
% uiwait(gcf);
if gN==7
    annularPar=eval(get(handles.pb_g7_et,'string'));
    y0=annularPar(1);
    x0=annularPar(2);
    r1=annularPar(3);
    r2=annularPar(4);
    r1_tmp=min([r1,r2]);
    r2_tmp=max([r1,r2]);
    r1=r1_tmp;
    r2=r2_tmp;
%     circumference=2*pi*r1;
%     pixel
    section=360;
    if 2*pi*r2>360
        section=round(2*pi*r2);
    end
    theta=linspace(0,2*pi,section);
    y1=sin(theta)*r1+y0;
    x1=cos(theta)*r1+x0;
    y2=sin(-theta)*r2+y0;
    x2=cos(-theta)*r2+x0;
    cx=[x1,x2,x1(end)].';
    cy=[y1,y2,y1(end)].';
else
    if length(pointPst)<gN
        pointPst{gN}=[];
%     handles.pointPst=cell(groupNO,1);
    handles.pointPstH{gN}=[];
    end
    if get(handles.cb_ROIuserInput,'value')
        dataPosition=get(handles.t_positions,'data');
        dataPosition=dataPosition(1:2,1:2);
        dataPosition2=dataPosition;
        dataTrash=zeros(2,2);
        dataTrash(1,1)=eval(dataPosition{1,1});
        dataTrash(1,2)=eval(dataPosition{1,2});
        dataTrash(2,1)=eval(dataPosition{2,1});
        dataTrash(2,2)=eval(dataPosition{2,2});
        dataPosition=dataTrash;
%         dataPosition=cell2mat(dataPosition);
        cy=[dataPosition(1,1)-.1,dataPosition(2,1)+.1,dataPosition(2,1)+.1,dataPosition(1,1)-.1]';
        cx=[dataPosition(1,2)-.1,dataPosition(1,2)-.1,dataPosition(2,2)+.1,dataPosition(2,2)+.1]';
    else
        if get(handles.cb_circle,'value')
            htrash=imellipse;
            
            addNewPositionCallback(htrash,@(p) title(mat2str(p,3)));
            fcn = makeConstrainToRectFcn('imellipse',get(gca,'XLim'),get(gca,'YLim'));
            setPositionConstraintFcn(htrash,fcn);
            %     wait( hAll{ii});
            %                 h = imellipse;
            vertices = wait(htrash);
            cx=vertices(:,1);
            cy=vertices(:,2);
            delete(htrash);
        else
%             [cx,cy,c]=improfile;
            [c,cx,cy]=roipoly;
        end
    if length(cx)==1
        cx=round(cx);
        cy=round(cy);
        cx=[cx-.1,cx-.1,cx+.1,cx+.1].';
        cy=[cy-.1,cy+.1,cy+.1,cy-.1].';
    end
    end
end
if get(handles.cb_circle2,'value')
    img=handles.img;
   bw=roipoly(img(:,:,1),cx,cy);
    t=regionprops(bw,'centroid');
    ax=t.Centroid(1);
    ay=t.Centroid(2);
    rr=(cx-ax).^2+(cy-ay).^2;
    r2=sqrt(rr);
    r=min(r2);
    r=r(1);
    cx=(ax+cos(linspace(0,2*pi,100))*r).';
    cy=(ay+sin(linspace(0,2*pi,100))*r).';

end
if get(handles.cb_gr11,'value')
    
    if length(handles.pointPst)>=1
        cx2=handles.pointPst{1};
        if length(cx2)>=1
            cx3=handles.pointPst{1}{1};
            if length(cx3)>=1
                img=handles.img;
               bw=roipoly(img(:,:,1),cx,cy);
                t=regionprops(bw,'centroid');
                ax=t.Centroid(1);
                ay=t.Centroid(2);
               bw2=roipoly(img(:,:,1),cx3(:,1),cx3(:,2));
                t2=regionprops(bw2,'centroid');
                ax2=t2.Centroid(1);
                ay2=t2.Centroid(2);
                
                cx=cx3(:,1)+ax-ax2;
                cy=cx3(:,2)+ay-ay2;
            end
        end
    end
end
data=pointPst{gN};
if isempty(data)
    data=cell(1,1);
    data{1}=[cx,cy];
else
    [m,n]=size(data);
    data{m+1,1}=[cx,cy];
end
pointPst{gN}=data;
handles.pointPst=pointPst;
axes(handles.axes_rawImg);
% set(gcf,'currentAxes',handles.axes_rawImg)
% zoom out;
hold on;
% if get(handles.cb_ROIuserInput,'value') && abs(max(cy)-min(cy))<1 && abs(max(cx)-min(cx))<1
if abs(max(cy)-min(cy))<1 && abs(max(cx)-min(cx))<1
%     get(handles.cb_ROIuserInput,'value') && abs(cy(2)-cy(1))<1 && abs(cx(2)-cx(1))<1
    hh=plot(cx(1),cy(1),'+','color',colorTmp(gN,:));
else
    hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(gN,:));
end
if isfield(handles, 'pointCurrent_h')~=1
handles.pointCurrent_h=hh;
handles.pointCurrent_group=gN;
else
    handles.pointCurrent_h=[handles.pointCurrent_h;hh];
    handles.pointCurrent_group=[handles.pointCurrent_group;gN];
end
if isfield(handles, 'pointPstH')~=1
handles.pointPstH=cell(size(handles.pointPst));
else
    handles.pointPstH{gN}=[handles.pointPstH{gN},hh];
end
handles=GroupROITag(handles);
function handles=GroupROITag(handles)
pointPst=handles.pointPst;
if isfield(handles,'grTag')
    grTag=handles.grTag;
    for ii=1:length(grTag)
        for jj=1:length(grTag{ii})
            if ishandle(grTag{ii}{jj})
                delete(grTag{ii}{jj})
            end
        end
        
    end
end
grTag={};
img=differentTypeRead(fullfile(handles.filePath,handles.file),handles.fileType);
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
%             vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            t=regionprops(bw,'centroid');
            ax=t.Centroid(1);
            ay=t.Centroid(2);
            grTag{ii}{jj}=text(ax,ay,['(',num2str(ii),',',num2str(jj),')'],'color',handles.colorTmp(ii,:));
%             ROItmp=img(bw);
%             rawData(kk,vv)=mean(ROItmp(:));
        end
    end
    handles.grTag=grTag;

function pushbutton12_Callback(hObject, eventdata, handles)
cd(handles.filePath)
[file,filePath]=uiputfile('parameterMat.mat');
cd(filePath)
%pointCurrent_h=handles.pointCurrent_h;
% pointCurrent_group=handles.pointCurrent_group;
pointPst=handles.pointPst;
colorTmp=handles.colorTmp;
save(file,'pointPst','colorTmp')
fileTmp=file;
filePathTmp=filePath;
filePath=handles.filePath;
fileType=handles.fileType;
cd(filePath)
img=differentTypeRead(handles.file,fileType);

[mm,nn,nnTrash]=size(img);
newImgColor=handles.newImgColor;
figure(78);imshow(newImgColor);hold on;
% handles=GroupROITag(handles);
 cd(filePathTmp)
grTag={};
for ii=1:length(pointPst)
    mapGroup=false(mm,nn);
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        plot(cx,cy,'color',handles.colorTmp(ii,:))
        
        bw=roipoly(img(:,:,1),cx,cy);
                   t=regionprops(bw,'centroid');
            ax=t.Centroid(1);
            ay=t.Centroid(2);
            grTag{ii}{jj}=text(ax,ay,['(',num2str(ii),',',num2str(jj),')'],'color',handles.colorTmp(ii,:));
        mapGroup=or(mapGroup, bw);
%         pointPst{ii}{jj}(:,1)=pointPst{ii}{jj}(:,1)+3;
    end
    imwrite(mapGroup,[fileTmp(1:end-4),'Group',num2str(ii),'.tif'],'tif','compression','none');
end
 saveas(78,'groupROITag.tif','tiff');
guidata(hObject,handles)
function pushbutton13_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.mat','MultiSelect','on');
cd(filePath)
if iscell(file)==0
    load(file)
else
    p=length(file);
    maxG=zeros(p,1);
    for ii=1:p
        load(file{ii});
        maxG(ii)=length(pointPst);
    end
    pointPstTmp=cell(p,1);
    for ii=1:p
        load(file{ii});
        for jj=1:maxG(ii)
            pointPstTmp{jj}=[pointPstTmp{jj}; pointPst{jj}];
        end
    end
    pointPst=pointPstTmp;
end
handles.pointPst=pointPst;
handles.colorTmp=colorTmp;
handles.pointCurrent_h=zeros(0,0);
handles.pointCurrent_group=zeros(0,0);
axes(handles.axes_rawImg)

% set(gcf,'currentAxes',handles.axes_rawImg)
hold on;
pointPstH=cell(size(pointPst));
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
        handles.pointCurrent_h=[handles.pointCurrent_h;hh];
        handles.pointCurrent_group=[handles.pointCurrent_group;ii];
        pointPstH{ii}=[pointPstH{ii},hh];
    end
end
handles.pointPstH=pointPstH;
figure;
cd(handles.filePath)
img=differentTypeRead(handles.file,handles.fileType);
mapRange=eval(get(handles.et_fileRange,'string'));
contents = cellstr(get(handles.pp_fileColormap,'String'));
colorSelected=contents(get(handles.pp_fileColormap,'value'));
imshow(img,mapRange);
colormap(eval(colorSelected{1}))
colorbar
hold on;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
%         handles.pointCurrent_h=[handles.pointCurrent_h;hh];
%         handles.pointCurrent_group=[handles.pointCurrent_group;ii];
    end
end
axes(handles.axes_rawImg)
handles=GroupROITag(handles);
guidata(hObject,handles)
% --- Executes on button press in normCorrwithinROI.
function normCorrwithinROI_Callback(hObject, eventdata, handles)
% hObject    handle to normCorrwithinROI (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
ROI_mean=cell(groupN,1);
pearsonMap=cell(groupN,1);
imageMap=zeros(m,n);
%%
if get(handles.et_meanFirstROI,'value')>=3 && get(handles.et_meanFirstROI,'value')<=8
    mm_th=get(handles.et_meanFirstROI,'value')-2;
     h_wait=waitbar(0,'waitbarName');
     referenceSignal=zeros(p,1,'single');
     for kk=1:p
         waitbar(kk/p,h_wait,['meanOfGroup',num2str(mm_th),'complete',num2str(100*kk/p,'%03.1f')]);
         file=fileName(kk).name;
         img=differentTypeRead(file,fileType);
         ii=mm_th;
         for jj=1:length(pointPst{mm_th})
             cx=pointPst{ii}{jj}(:,1);
             cy=pointPst{ii}{jj}(:,2);
             bw=roipoly(img,cx,cy);
             referenceSignal(kk)=referenceSignal(kk)+sum(img(bw));
         end
     end
     mm_thTotal=0;
     ii=mm_th;
     for jj=1:length(pointPst{mm_th})
         cx=pointPst{ii}{jj}(:,1);
         cy=pointPst{ii}{jj}(:,2);
          bw=roipoly(img,cx,cy);
          mm_thTotal=mm_thTotal+sum(bw(:));
     end
     referenceSignal=referenceSignal/mm_thTotal;
     close(h_wait)
     handles.referenceSignal=referenceSignal;
end
%%
for ii=1:length(pointPst)
    ROI_mean{ii}=zeros(p,ROI_N(ii),'single');
    pearsonMap{ii}=cell(ROI_N(ii),1);
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        bw=roipoly(img,cx,cy);
        n=sum(bw(:));
        ROI=zeros(p,n,'single');
        % ROI=zeros(p,n);
        cd(filePath)
        waitbarName=['group',num2str(ii),'ROI',num2str(jj),'complete'];
        h_wait=waitbar(0,'waitbarName');
        for kk=1:p
            waitbar(kk/p,h_wait,[waitbarName,num2str(100*kk/p,'%03.1f')]);
            file=fileName(kk).name;
            img=differentTypeRead(file,fileType);
%             bw=roipoly(img,cx,cy);
            ROI(kk,:)=img(bw);
        end
         ROI=IOS_time_gui_filter(ROI,handles,'onYdirection');
        close(h_wait)
        ROI_mean{ii}(:,jj)=mean(ROI,2);
        %%
        %% get pearson correlation coefficients
        pearsonMap{ii}{jj}=zeros(n,1,'single');
        if get(handles.et_meanFirstROI,'value')==2 % mean of first ROI as reference
            for ss=1:n
               pearsonMap{ii}{jj}(ss)=corr(ROI_mean{1}(:,1),ROI(:,ss));
               handles.referenceSignal=ROI_mean{1}(:,1);
            end
        elseif get(handles.et_meanFirstROI,'value')==1 % mean of each ROI as reference
            for ss=1:n
               pearsonMap{ii}{jj}(ss)=corr(ROI_mean{ii}(:,jj),ROI(:,ss));
            end
        elseif get(handles.et_meanFirstROI,'value')>=3 && get(handles.et_meanFirstROI,'value')<=8% mean of mm_th group as reference
           mm_th=get(handles.et_meanFirstROI,'value')-2;
            for ss=1:n
               pearsonMap{ii}{jj}(ss)=corr(referenceSignal,ROI(:,ss));
            end
        end
        %% map correlation coefficient back to image
        imageMap(bw)=pearsonMap{ii}{jj};
    end
end
handles.imageMap=imageMap;
handles.ROI_mean=ROI_mean;
% mapRange=eval(get(handles.et_caxis,'string'));
% contents = cellstr(get(handles.pp_colormap,'String'));
% colorSelected=contents(get(handles.pp_colormap,'value'));
% figure; imagesc(imageMap,mapRange);
% colormap(eval(colorSelected))
guidata(hObject,handles)
function et_pearsonFolder_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
function et_pearsonFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_pearsonFolder (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in cb_inverseNegative.
function cb_inverseNegative_Callback(hObject, eventdata, handles)
% --- Executes on selection change in pp_colormap.
function pp_colormap_Callback(hObject, eventdata, handles)
handles=pearsonMapShow(handles);
guidata(hObject,handles)
function pp_colormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_caxis_Callback(hObject, eventdata, handles)
% hObject    handle to et_caxis (see GCBO)
handles=pearsonMapShow(handles);
guidata(hObject,handles)
function et_caxis_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton15_Callback(hObject, eventdata, handles)
handles=pearsonMapShow(handles);
guidata(hObject,handles)
function handles=pearsonMapShow(handles)
if isfield(handles,'imageMap')
imageMap=handles.imageMap;
mapRange=eval(get(handles.et_caxis,'string'));
contents = cellstr(get(handles.pp_colormap,'String'));
colorSelected=contents(get(handles.pp_colormap,'value'));
figure(99); imshow(imageMap,mapRange);
colormap(eval(colorSelected{1}));
handles.colorSelectedRange=eval(colorSelected{1});
colorbar;
hold on;
if isfield(handles,'colorTmp') && isempty(handles.pointPst)==0
    colorTmp=handles.colorTmp;
    pointPst=handles.pointPst;
    for ii=1:1%ii=1:length(pointPst)
        for jj=1:1%jj=1:length(pointPst{ii})
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(3,:));
%             hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
%             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
%             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
        end
    end
end
hold off;
%% correlation cross groups
ROI_mean=handles.ROI_mean;
p1=length(ROI_mean);
p2=length(ROI_mean{1}(:,1));
data=zeros(p2,p1);
colorTmp=handles.colorTmp;
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p2).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
figure(100);
for ii=1:p1
    if ii>=2
        hold on;
    end
    data(:,ii)=mean(ROI_mean{ii},2);
    plot(timeCourse,data(:,ii),'color',colorTmp(ii,:))
   % plot(linspace(10.1,60,length(data)),data(:,ii),'color',colorTmp(ii,:))
end
global ttt
ttt=data(:,1);
grid on;
hold off;
xlabel(timeCourseUnit)
disp('pairwise Pearson correlation between groups%%%%%%%%%%%');
global correlationBtnGroups
correlationBtnGroups=corr(data)
handles.correlationBtnGroups=correlationBtnGroups;
handles.meanOfGroup=data;
end
function pushbutton16_Callback(hObject, eventdata, handles)
% --- Executes on button press in pb_IOScolor.
function pb_IOScolor_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
%%
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1=get(handles.et_rawIOS_gray,'string');
resultName2=get(handles.et_rawIOS_pseudocolor,'string');
resultPath1=fullfile(filePathTmp,resultName1);
resultPath2=fullfile(filePathTmp,resultName2);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
if exist(resultPath2)==7
else
    mkdir(resultName2);
end
cd(filePath)
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
fileName=handles.fileName;
p=length(fileName);
newIntensity=inline(get(handles.et_IOSinline,'string'));
contents = cellstr(get(handles.pp_IOScolormap,'String'));
colorSelectedTmp=contents(get(handles.pp_IOScolormap,'value'));
colorSelected=eval(colorSelectedTmp{1});
hh=waitbar(0,'please wait ...');
uint16Flag=get(handles.cb_unit16,'value');
for ii=1:p
                waitbar(ii/p,hh, ...
                 [num2str(100*ii/p,'%03.1f'),'% completed']);
    file=fileName(ii).name;
    cd(filePath)
    img=differentTypeRead(file,fileType);
    img2=newIntensity(img);
    if uint16Flag~=1
    img3=uint8(img2);
    elseif uint16Flag==1
      img3=uint16(img2);
    end
    cd(resultPath1)
    saveName1=[file(1:end-length(fileType)),'tif'];
    imwrite(img3,saveName1,'tiff','compression','none');
    cd(resultPath2)
    saveName2=[file(1:end-length(fileType)),'png'];
    imgColor=ind2rgb(img3,colorSelected);
    imwrite(imgColor,saveName2,'png')
    if get(handles.cb_movie,'value')
        M_gray(ii)=im2frame(img3,gray(256));
        M_color(ii)=im2frame(img3,colorSelected);
    end
end
 if get(handles.cb_movie,'value')
    movie2avi(M_gray,fullfile(resultPath2,'movie'),'compression','none','fps',5);
    implay(M_gray)
    movie2avi(M_color,fullfile(resultPath2,'movie'),'compression','none','fps',5);
    % movie2avi(M2,fullfile(filePath,'movie2'),'compression','none','fps',5);
    implay(M_color)
 end
close(hh)
guidata(hObject,handles)
function et_rawIOS_gray_Callback(hObject, eventdata, handles)
function et_rawIOS_gray_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_rawIOS_pseudocolor_Callback(hObject, eventdata, handles)
function et_rawIOS_pseudocolor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_IOSinline_Callback(hObject, eventdata, handles)
function et_IOSinline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_IOSinline (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pp_IOScolormap.
function pp_IOScolormap_Callback(hObject, eventdata, handles)
function pp_IOScolormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_binning_Callback(hObject, eventdata, handles)
function et_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton19_Callback(hObject, eventdata, handles)
cd(handles.filePath)
result=get(handles.et_pearsonFolder,'string');
resultFolder=fullfile(handles.filePath,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
[file,filePath]=uiputfile('imageRaw.mat');
cd(filePath)
%pointCurrent_h=handles.pointCurrent_h;
% pointCurrent_group=handles.pointCurrent_group;
imageMap=handles.imageMap;
mapRange=eval(get(handles.et_caxis,'string'));
contents = cellstr(get(handles.pp_colormap,'String'));
colorSelected=contents(get(handles.pp_colormap,'value'));
a=mapRange(1);
b=mapRange(2);
imageMapSave=(b-256)/(a-1)*(imageMap-a)+1;
imgColor=ind2rgb(uint8(imageMapSave),eval(colorSelected{1}));
imwrite(imgColor,'correlationMap_nonCompression.tif','tif','compression','none')
ROI_mean=handles.ROI_mean;
colorTmp=handles.colorTmp;
correlationBtnGroups=handles.correlationBtnGroups;
save(file,'imageMap','colorTmp','ROI_mean')
save('correlationBtnGroups.txt','correlationBtnGroups','-ASCII')
figure(99)
saveas(gcf,'corrWithinGroup','tiff')
figure(100)
saveas(gcf,'meanOfGroups','tiff')
meanOfGroup=handles.meanOfGroup;
[m,n]=size(meanOfGroup);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),m).';
meanOfGroupSaver=zeros(m,n+1);
meanOfGroupSaver(:,1)=timeCourse;
meanOfGroupSaver(:,2:end)=meanOfGroup;
save('meanOfGroups.txt','meanOfGroupSaver','-ASCII')
guidata(hObject,handles)
function pushbutton20_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.mat');
cd(filePath)
load(file)
handles.imageMap=imageMap;
handles.ROI_mean=ROI_mean;
handles.colorTmp=colorTmp;
guidata(hObject,handles)
% --- Executes on selection change in et_meanFirstROI.
function et_meanFirstROI_Callback(hObject, eventdata, handles)
function et_meanFirstROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_meanFirstROI (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pp_fileColormap.
function pp_fileColormap_Callback(hObject, eventdata, handles)
handles=singleColorImageShow(handles);
guidata(hObject,handles)
function pp_fileColormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_fileColormap (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_fileRange_Callback(hObject, eventdata, handles)
% hObject    handle to et_fileRange (see GCBO)
handles=singleColorImageShow(handles);
guidata(hObject,handles)
% et_fileRange as text
%        str2double(get(hObject,'String')) returns contents of et_fileRange as a double
function et_fileRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_fileRange (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton21_Callback(hObject, eventdata, handles)
handles=singleColorImageShow(handles);
guidata(hObject,handles)
function handles=singleColorImageShow(handles)
figure(98);
cd(handles.filePath)
img=differentTypeRead(handles.file,handles.fileType);
mapRange=eval(get(handles.et_fileRange,'string'));
contents = cellstr(get(handles.pp_fileColormap,'String'));
colorSelected=contents(get(handles.pp_fileColormap,'value'));
imshow(img,mapRange);
colormap(eval(colorSelected{1}))
ht=colorbar;
hold on;
if isfield(handles,'colorTmp') && isempty(handles.pointPst)==0
    colorTmp=handles.colorTmp;
    pointPst=handles.pointPst;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            if abs( abs(cy(2)-cy(1))-.2)<.00001 && abs(cx(2)-cx(1))==0
                hh=plot(cx(1),cy(1),'+','color',colorTmp(ii,:));
            else
%             handles.pointPst{ii}{jj}(:,1)=pointPst{ii}{jj}(:,1)-14;
%             handles.pointPst{ii}{jj}(:,2)=pointPst{ii}{jj}(:,2)+9;
            hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
            end
%             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
%             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
        end
    end
end
hold off;
function dt_IntestedImg_Callback(hObject, eventdata, handles)
% hObject    handle to dt_IntestedImg (see GCBO)
cd(handles.filePath)
handles.fileName=dir(['*.',handles.fileType]);for iss=length(handles.fileName):-1:1;if strcmp(handles.fileName(iss).name(1:2),'._'); handles.fileName(iss)=[];end;end
interestedImg=get(hObject,'String');
tt=eval(interestedImg);
handles.fileName=handles.fileName(tt);
handles.p=length(handles.fileName);
handles.file=handles.fileName(1).name;
imageShow(handles);
%show images' name
p=handles.p;
fileName=handles.fileName;
listName=cell(p,1);
for ii=1:p
    listName{ii}=fileName(ii).name;
end

handles.oriRange = [tt(1),tt(end)];
set(handles.lb_curFileNames,'value',tt(1));
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
guidata(hObject,handles);
% dt_IntestedImg as text
%        str2double(get(hObject,'String')) returns contents of dt_IntestedImg as a double
function dt_IntestedImg_CreateFcn(hObject, eventdata, handles);
% hObject    handle to dt_IntestedImg (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton22_Callback(hObject, eventdata, handles)
% --- Executes on button press in cb_ROIuserInput.
function cb_ROIuserInput_Callback(hObject, eventdata, handles)
% hObject    handle to cb_ROIuserInput (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of cb_ROIuserInput
% --- Executes on button press in cb_movie.
function cb_movie_Callback(hObject, eventdata, handles)
% hObject    handle to cb_movie (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of cb_movie
function pushbutton24_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
% cd(filePath)
img=differentTypeRead(fullfile(filePath,file),fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
rawData=zeros(p,sum(ROI_N(:))+1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
h_wait=waitbar(0,'wait');
for kk=1:p
    %h_wait=waitbar(0,'wait');
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
 if get(handles.pp_spfilterMethod,'value')>1
    img=differentTypeReadFilter_handles(img,handles);
 end     
    vv=1;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
            vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            ROItmp=img(bw);
            rawData(kk,vv)=mean(ROItmp(:));
        end
    end
end
close(h_wait)
rawData(:,2:end)= ...
    IOS_time_gui_filter(rawData(:,2:end),handles,'vertical');
% rawData=log(rawData);
handles.rawData_meanOfROI=rawData;
guidata(hObject,handles)
function pushbutton25_Callback(hObject, eventdata, handles)
pathDefault=cd;
% cd(handles.filePath)
result=['meanROI'];
% if get(handles.pp_spfilterMethod,'value')>1
%     result=[result,handles.filterMethod.saveName];
% end
%  
resultFolder=fullfile(handles.filePath,result);
% if exist(resultFolder)==7
% else
    mkdir(handles.filePath,result)
% end
cd(resultFolder)
[file,filePath]=uiputfile('meanOfROI.txt');
cd(filePath)
meanOfROI=handles.rawData_meanOfROI;
figure(102)
saveas(gcf,'meanofROI','tiff')
save(file,'meanOfROI','-ASCII');
%halfT=handles.halfT;
%maxV=halfT.maxV;
%peakTime=halfT.peakTime;
%halfTime=halfT.halfTime;
file2=file(1:end-4);
%save([file,'maxV.txt'],'maxV','-ASCII');
%save([file,'peakTime.txt'],'peakTime','-ASCII');
%save([file,'halfTime.txt'],'halfTime','-ASCII');
% [file,filePath]=uiputfile('parameterMat.mat');
% cd(filePath)
%pointCurrent_h=handles.pointCurrent_h;
% pointCurrent_group=handles.pointCurrent_group;
pointPst=handles.pointPst;
colorTmp=handles.colorTmp;
save([file,'_ROI_Parameters.mat'],'pointPst','colorTmp')
%meanOfROI_IOS=handles.rawData_meanOfROI_IOS;
% cd(pathDefault)
%figure(144)
%file='meanOfROI_IOS.txt';
%saveas(gcf,'meanofROI_IOS','tiff')
%save(file,'meanOfROI_IOS','-ASCII');
function et_timeCourse_Callback(hObject, eventdata, handles)
% hObject    handle to et_timeCourse (see GCBO)
% et_timeCourse as text
%        str2double(get(hObject,'String')) returns contents of et_timeCourse as a double
function et_timeCourse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_timeCourse (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_unit_Callback(hObject, eventdata, handles)
% hObject    handle to et_unit (see GCBO)
% et_unit as text
%        str2double(get(hObject,'String')) returns contents of et_unit as a double
function et_unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_unit (see GCBO)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton26_Callback(hObject, eventdata, handles)
rawData=handles.rawData_meanOfROI;
rawDatadIOS=handles.rawData_meanOfROI;
rawDatadI=handles.rawData_meanOfROI;

p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(102);hold off;
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-*','--',':','-.','o-','*--','*:','*-.'};

% IOS
pre=eval(get(handles.et_pre2,'string'));
if pre>p
    pre=1;
end

rawDatadIOS(:,2:end)=dIOSization(rawData(:,2:end),pre);
rawDatadI(:,2:end)=dIization(rawData(:,2:end),pre);

vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        subplot(3,1,1);
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'LineWidth',2,'MarkerSize',2,'color',colorTmp(ii,:))
        title('F','fontsize',15)
        set(gca, 'FontSize', 15);
        grid on;
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
       
              ccc=(rawData(:,vv));snrCCC=mean(ccc)/std(ccc);
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'SNR:',num2str(snrCCC)];

    end
end
hold off

vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        subplot(3,1,2);
        plot(rawDatadIOS(:,1),rawDatadIOS(:,vv),ROIii{jj},'LineWidth',2,'MarkerSize',2,'color',colorTmp(ii,:))
        title('\DeltaF/F','fontsize',15)
        set(gca, 'FontSize', 15);
        grid on;
    end
end
hold off

vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        subplot(3,1,3);
        plot(rawDatadI(:,1),rawDatadI(:,vv),ROIii{jj},'LineWidth',2,'MarkerSize',2,'color',colorTmp(ii,:))
        title('\DeltaF','fontsize',15)
        set(gca, 'FontSize', 15);
    end
end




h_legend=legend(legendName);
set(h_legend,'FontSize',16);
% set(gca,'fontsize',15)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;
handles.rawData_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;

%% IOS
% pre=eval(get(handles.et_pre2,'string'));
% if pre>p
%     pre=1;
% end
% rawData(:,2:end)=dIOSization(rawData(:,2:end),pre);
% %
% [halfT.maxV,halfT.peakTime,halfT.halfTime]=halfTime_cal(handles.rawData_meanOfROI,pre);
% handles.halfT=halfT;
% %
% legendName=cell(size(rawData,2)-1,1);
% figure(144);hold off;
% vv=1;
% for ii=1:length(pointPst)
%     for jj=1:length(pointPst{ii})
%         vv=vv+1;
%         if vv>2
%             hold on;
%         end
%         plot(rawData(:,1),rawData(:,vv),ROIii{jj},'LineWidth',4,'MarkerSize',4,'color',colorTmp(ii,:))
%        %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
% %               ccc=(rawData(:,vv));snrCCC=mean(ccc)/std(ccc);
% %         legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'halfT',num2str(halfT.halfTime(vv-1),'%6.3f'), ...
% %             'peakT',num2str(halfT.peakTime(vv-1),'%6.3f'), ...
% %             'maxV',num2str(halfT.maxV(vv-1),'%6.3f'),];
%         
%     end
% end
% h_legend=legend(legendName);
% set(h_legend,'FontSize',12);
% set(gca,'fontsize',15)
% title('\DeltaI/I','fontsize',15)
% % xLimit=get(gca,'xlim');
% % yLimit=get(gca,'ylim');
% % pre2=eval(get(handles.et_pre2,'string'));
% % line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% % set(gca,'xlim',xLimit)
% % set(gca,'ylim',yLimit)
% grid on;
% handles.rawData_meanOfROI_IOS=rawData;


%% obtian 
guidata(hObject,handles)
function [maxV,peakTime,halfTime]=halfTime_cal(varargin)
% this function calculate the max value, time to peak and time to half
% peak;
% try: [maxV,peakTime,halfTime]=halfTime(t,data,pre);
% if the syntax [maxV,peakTime,halfTime]=halfTime(data,pre) is used, then the
% first column of data is taken as time; if data only has one column, then
% t equals to 1:length(data);
%if the syntax [maxV,peakTime,halfTime]=halfTime(data) is used, then pre=1;
%if there is no input argument, then the user is offered to select the
%data;
% check input aruments
if nargin<=2
    if nargin==0
        [file,filePath]=uigetfile('*.*');
        data=load(fullfile(filePath,file));
        pre=1;
    elseif nargin==1
        data=varargin{1};
        pre=1;
    elseif nargin==2
        data=varargin{1};
        pre=varargin{2};
    end
    [m,~]=size(data);   
    if m==1
        data=data.';
        
    end
    [m,n]=size(data); 
    if n==1
        t=(1:m).';
    else
        t=data(:,1);
    end
    data2=data(:,2:end);

    

    elseif nargin==3
        t=varargin{1}(:);
        data2=varargin{2}; 
        pre=varargin{3};  
        
else
    return;
    
end

% IOS?
data2=dIOSization(data2,pre(end));

p0=size(data2,2);
peakTime=zeros(p0,1);
halfTime=zeros(p0,1);
maxV=zeros(p0,1);
halfPeak=zeros(p0,1);
figure(3);hold off;
step=0;
for ii=1:p0
    
%     hold off;
    
    y=data2(:,ii);
    
    hold on;
    [Y,I]=max(y);
    I=I(1);
    Y=Y(1);
    maxV(ii)=Y;
    peakTime(ii)=t(I);
    % find halfTime
    halfPeak(ii)=(Y+mean(y(pre)))/2;
    y1=y(1:I);
    y2=abs(y1-halfPeak(ii));
    [Y2,I2]=min(y2);
    Y2=Y2(1);
    I2=I2(1);
    if Y2==0
        halfTime(ii)=t(I2);
    else
        nextMin=y2(I2);
        if nextMin<halfPeak(ii)
            t1=t(I2); t2=t(I2+1);
            x1=y1(I2);x2=y1(I2+1);
%         else
%             t1=t(I2-1); t2=t(I2);
%             x1=y1(I2-1);x2=y1(I2);
        end
        P = polyfit([x1,x2],[t1,t2],1);
        halfTime(ii)=polyval(P,halfPeak(ii));
        
    end

end
step=0;
for ii=1:p0
    step2=max(maxV(:))-mean(y(1:pre));
    y=data2(:,ii);
    plot(t,y-step,'-.');
    line([t(1),t(end)],halfPeak(ii)*[1,1]-step);
    
    plot(halfTime(ii),halfPeak(ii)-step,'o')
    plot(peakTime(ii),maxV(ii)-step,'*')
    xlabel('t','fontsize',15)
    ylabel('Intensity','fontsize',15)
    step=step2+step;    
end

figure(4);
subplot(1,3,1);
plot(1:p0,maxV,'r-o');
xlabel('ROI','fontsize',15)
ylabel('maxV','fontsize',15)

subplot(1,3,2);
plot(1:p0,peakTime,'r-o');
xlabel('ROI','fontsize',15)
ylabel('peakTime','fontsize',15)
subplot(1,3,3);
plot(1:p0,halfTime,'r-o');
xlabel('ROI','fontsize',15)
ylabel('halfTime','fontsize',15)
function pushbutton27_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.txt');
cd(filePath)
handles.rawData_meanOfROI=load(file);
guidata(hObject,handles)
function pushbutton30_Callback(hObject, eventdata, handles)
figure(eval(get(handles.et_fig,'string')))
xLimit=get(gca,'xlim');
yLimit=get(gca,'ylim');
pre2=eval(get(handles.et_verticalLine,'string'));
hold on;
line([pre2,pre2],yLimit,'color','k')
set(gca,'xlim',xLimit)
set(gca,'ylim',yLimit)
hold off;
function et_verticalLine_Callback(hObject, eventdata, handles)
% hObject    handle to et_verticalLine (see GCBO)
% et_verticalLine as text
%        str2double(get(hObject,'String')) returns contents of et_verticalLine as a double
function et_verticalLine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_verticalLine (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_fig_Callback(hObject, eventdata, handles)
function et_fig_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pp_tpfilterMethod.
function pp_tpfilterMethod_Callback(hObject, eventdata, handles)
handles=pp_filterMethod_outside(hObject,handles);
guidata(hObject,handles);
function handles=pp_filterMethod_outside(hObject,handles)
handles.filterMethod.name=get(hObject,'Value')-1;
handles=filterInfoRead2(handles);
nameTmp=handles.filterMethod.name;
if nameTmp==0 %non filter
    set(handles.p_nonFilter,'visible','off');
elseif nameTmp==1% gaussian filter
    set(handles.p_nonFilter,'visible','on');
    set(handles.p_wiener,'visible','on')
elseif nameTmp==2 % median filter
   set(handles.p_nonFilter,'visible','on');
   set(handles.p_wiener,'visible','off');
elseif nameTmp==3
    set(handles.p_nonFilter,'visible','on');
    set(handles.p_wiener,'visible','off');
elseif nameTmp==4
    set(handles.p_nonFilter,'visible','off');
    prompt = {'fitting order:','fitting section:','baseline(0 or mean):'};
    dlg_title = 'Input for detrend function';
    num_lines = 1;
    def = {'2',['1:',num2str(handles.p)],'mean'};
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    handles.filterMethod.detrend_order=eval(answer{1});
    handles.filterMethod.detrend_section=eval(answer{2});
    handles.filterMethod.detrend_base=answer{3};
end
% guidata(hObject,handles)
function handles=filterInfoRead2(handles)
handles.filterMethod.sizNo=str2num(get(handles.et_sizNo,'string'));
%handles.sizSuper=str2num(get(handles.et_sizSuper,'string'));
handles.filterMethod.sigmaNo=str2num(get(handles.et_sigmaNo,'string'));
%handles.sigmaSuper=str2num(get(handles.et_sigmaSuper,'string'));
handles.filterMethod.name=get(handles.pp_tpfilterMethod,'Value')-1;

function t1=IOS_time_gui_filter(t1,handles,varargin)
handles=filterInfoRead2(handles);
filterMethod=handles.filterMethod;
sizNo=filterMethod.sizNo;
sigmaNo=filterMethod.sigmaNo;
nameNo=filterMethod.name;
siz=sizNo;
sigma=sigmaNo;
width=siz;
[m1,n1]=size(t1);
if nargin==2
    if n1>m1
        t1=t1.';
    end
    if nameNo==0
        %imgAvg=imgAvg;
    elseif nameNo==1
        w=fspecial('gaussian',[siz,1],sigma);
        t1=imfilter(t1,w,'replicate');
    elseif nameNo==2
        t1=IOS_time_medianFilterAdapt(t1,siz);
    elseif nameNo==3
            t1=wiener2(t1,[siz,1]);
    elseif nameNo==4
            detrend_order=handles.filterMethod.detrend_order;
            detrend_section=handles.filterMethod.detrend_section;
            detrend_base=handles.filterMethod.detrend_base;
            [m2,n2]=size(t1);
            if strcmp(detrend_base,'mean')        
                for ii=1:n2
                    t1(:,ii)=detrendnonlin(t1(:,ii),detrend_order,detrend_section);

                end
             else
                for ii=1:n2
                    t1(:,ii)=detrendnonlin(t1(:,ii),detrend_order,detrend_section,eval(detrend_base));

                end
            end
    end
    if n1>m1
        t1=t1.';
    end
elseif nargin==3
    if nameNo==0
        %imgAvg=imgAvg;
    elseif nameNo==1
        w=fspecial('gaussian',[siz,1],sigma);
        t1=imfilter(t1,w,'replicate');
    elseif nameNo==2
        t1=IOS_time_medianFilterAdapt(t1,siz);
    elseif nameNo==3
            t1=wiener2(t1,[siz,1]);
    elseif nameNo==4
            detrend_order=handles.filterMethod.detrend_order;
            detrend_section=handles.filterMethod.detrend_section;
            detrend_base=handles.filterMethod.detrend_base;
            [m2,n2]=size(t1);
            if strcmp(detrend_base,'mean')
                
            
                    for ii=1:n2
                        t1(:,ii)=detrendnonlin(t1(:,ii),detrend_order,detrend_section);
                        
                    end
 
                
            else
  
                    for ii=1:n2
                        t1(:,ii)=detrendnonlin(t1(:,ii),detrend_order,detrend_section,eval(detrend_base));
                        
                    end
    
                
            end
        
%         t1=edge(t1,'canny',[.05,.4],sigma);
    end
end
function x=IOS_time_medianFilterAdapt(x,n)
[m,n2]=size(x);
y=zeros(m+2*n,n2);
% y=zeros(m+2*n,n2,'single');
for ii=1:n2
    y(n+1:n+m,ii)=x(:,ii);
    y(1:n,ii)=x(1,ii);
    y(n+m+1:end,ii)=x(end,ii);
    y(:,ii)=medfilt1(y(:,ii),n);
    x(:,ii)=y(n+1:n+m,ii);
end
function pp_tpfilterMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_tpfilterMethod (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function b = detrendnonlin(a, varargin)
%% try detrendnonlin(a,order,section,baseline)
x=a;
%% input regulation
if isempty(varargin)
    order=5;
    section=1:length(x);
    baseline=mean(a(:));
elseif length(varargin)==1
    order=varargin{1};
    section=1:length(x);
    baseline=mean(a(:));  
elseif length(varargin)==2
    order=varargin{1};
    section=varargin{2};
    baseline=mean(a(:));  
elseif length(varargin)==3
    order=varargin{1};
    section=varargin{2};
    baseline=varargin{3};     
end


%% datafitting
a=a(:);
x=a(section);

p = polyfit((1:numel(x))', x(:), order);
yTemp=polyval(p, (1:numel(x))');
y = x(:) - yTemp;

b1=a(1:section(1)-1);
b2=y;
b3=a(section(end)+1:end);
if ~isempty(b1)
    b1=b1-b1(end)+b2(1)-(a(section(1))-a(section(1)-1));
end
if ~isempty(b3)
    b3=b3-b3(1)+b2(end)+(a(section(end)+1)-a(section(end)));
end
b=[b1;b2;b3]+baseline;

yTemp=yTemp;
figure(14);hold off; plot(1:length(a),a,'-b',section,yTemp,'k--',1:length(b),b,'r'); legend({'raw',['order=',num2str(order)],'detrended'})
grid on;
% hold on; plot(1:length(b),b,'g')
% grid on;
function et_sizNo_Callback(hObject, eventdata, handles)
function et_sizNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_sigmaNo_Callback(hObject, eventdata, handles)
function et_sigmaNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_NO_Callback(hObject, eventdata, handles)
function et_NO_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function handles=pushbutton31_Callback(hObject, eventdata, handles)
handles=IOS_time_parameterRead(handles);
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
fileType=handles.fileType;
pre1=handles.pre1;
pre2=handles.pre2;
preNO=pre2-pre1+1;
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
% [uV sV] = memory;
% ts=zeros(m,n,'uint16');
% ts_para=whos('ts');
% ts2=single(ts);
% ts2_para=whos('ts2');
if strcmp(fileType,'mat')~=1 %&& ts_para.bytes*preNO<uV.MaxPossibleArrayBytes*.9
    preData=zeros(m,n,preNO,'uint16');
    hh=waitbar(0,'wait...');
    for ii=pre1:pre2
            waitbar(ii/preNO,hh, ...
             [num2str(100*ii/preNO,'%03.1f'),'% completed']);
        img2=imread(fileName(ii).name);
        %img=img2(:,:,1);
        preData(:,:,ii)=(img2(:,:,1));
    end
    clear img2 img
    % for ii=1:m
    %     imgTmp=preData(ii,:,:);
    %
    % end
    close(hh)
    handles.meanImg=zeros(m,n,'single');
    handles.stdImg=zeros(m,n,'single');
    for ii=1:m
%         preDataTmp=squeeze(single(preData(ii,:,:)));
%         preDataTmp=IOS_time_gui_filter(preDataTmp,handles,'vertical');
        handles.meanImg(ii,:)=mean(single(preData(ii,:,:)),3);
        handles.stdImg(ii,:)=std(single(preData(ii,:,:)),0,3);
    end
    %%
    clear preData
elseif strcmp(fileType,'mat')==1 %&& ts2_para.bytes*preNO<uV.MaxPossibleArrayBytes*.9
  preData=zeros(m,n,preNO,'single');
    hh=waitbar(0,'wait...');
    for ii=pre1:pre2
            waitbar(ii/preNO,hh, ...
             [num2str(100*ii/preNO,'%03.1f'),'% completed']);
        img2=differentTypeRead(fileName(ii).name,fileType);
        %img=img2(:,:,1);
        preData(:,:,ii)=(img2(:,:,1));
    end
    clear img2 img
    % for ii=1:m
    %     imgTmp=preData(ii,:,:);
    %
    % end
    close(hh)
    handles.meanImg=zeros(m,n,'single');
    handles.stdImg=zeros(m,n,'single');
    for ii=1:m
%         preDataTmp=squeeze(single(preData(ii,:,:)));
%         preDataTmp=IOS_time_gui_filter(preDataTmp,handles,'vertical');
        handles.meanImg(ii,:)=mean(single(preData(ii,:,:)),3);
        handles.stdImg(ii,:)=std(single(preData(ii,:,:)),0,3);
    end
    %%
    clear preData
else
       handles.meanImg=zeros(m,n,'single');
    handles.stdImg=zeros(m,n,'single');
    hh=waitbar(0,'wait...');
    for ii=1:m
                    waitbar(ii/m,hh, ...
             [num2str(100*ii/m,'%03.1f'),'% completed']);
         preData=zeros(preNO,n,'single');
        for jj=pre1:pre2
            img=differentTypeRead(fileName(jj).name,fileType);
            preData(jj,:)=single(img(ii,:,1));
        end
        handles.meanImg(ii,:)=mean(preData);
        handles.stdImg(ii,:)=std(preData);
    end
    close(hh);
end
guidata(hObject,handles)

function handles=pushbutton32_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
fileName=handles.fileName;
method=get(handles.pm_method,'value');
if method==3 && isfield(handles,'meanImg')==0
    cd(filePath)
    img=differentTypeRead(file,fileType);
    handles.meanImg=img(:,:,1);
    handles.stdImg=img(:,:,1);
end
meanImg=handles.meanImg;
stdImg=handles.stdImg;
x1=eval(get(handles.et_x1,'string'));
x2=eval(get(handles.et_x2,'string'));
x3a=eval(get(handles.et_x3a,'string'));
x3b=eval(get(handles.et_x3b,'string'));
NO=round(eval(get(handles.et_NO,'string')));
p=length(fileName);
cd(filePath);
img=differentTypeRead(file,fileType);
[mm,nn,nn_trash]=size(img);
%% map
positiveMap=false(mm,nn);
negativeMap=false(mm,nn);
%%
% [uV sV] = memory;
% ts=zeros(mm,nn,'uint16');
% ts_para=whos('ts');
% zones=floor(ts_para.bytes*p/(uV.MaxPossibleArrayBytes*.9))+1;
pointPst=handles.pointPst;
cd(filePath)
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
 positiveNO=zeros(p,sum(ROI_N(:)),'single');
 negativeNO=zeros(p,sum(ROI_N(:)),'single');
positiveSum=zeros(p,sum(ROI_N(:)),'single');
negativeSum=zeros(p,sum(ROI_N(:)),'single');
 %positiveNegativeNO=zeros(p,sum(ROI_N(:)),'single');
%% to save time, doing calculation case by case
if get(handles.cb_ROIsmall,'value')==1 %% user define
 %%
 tk_ii=0;
 for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        bw=roipoly(img,cx,cy);
        ROImean=meanImg(bw);
        ROIstd=stdImg(bw);
        n=sum(bw(:));
        ROI=zeros(p,n,'single');
        % ROI=zeros(p,n);
        cd(filePath)
        waitbarName=['group',num2str(ii),'ROI',num2str(jj),'complete'];
        h_wait=waitbar(0,'waitbarName');
        for kk=1:p
            waitbar(kk/p,h_wait,[waitbarName,num2str(100*kk/p,'%03.1f')]);
            file=fileName(kk).name;
            img=differentTypeRead(file,fileType);
            %bw=roipoly(img,cx,cy);
            ROI(kk,:)=img(bw);
        end
         ROI=IOS_time_gui_filter(ROI,handles,'onYdirection');
        close(h_wait)
        tk_ii=tk_ii+1;
        [positiveNO(:,tk_ii),negativeNO(:,tk_ii), ...
            positiveMapROI,negativeMapROI,positiveSum(:,tk_ii),negativeSum(:,tk_ii)]= ...
            positiveNegative(ROI,ROImean.',ROIstd.',handles);
            I_positive=find(bw);
            bw_positive=false(mm*nn,1);
            bw_positive(I_positive(positiveMapROI))=1;
            positiveMap=or(positiveMap,reshape(bw_positive,size(positiveMap)));
            I_negative=find(bw);
            bw_negative=false(mm*nn,1);
            bw_negative(I_negative(negativeMapROI))=1;
            negativeMap=or(negativeMap,reshape(bw_negative,size(negativeMap)));
%         sum(bw(:))
% I_positive=find(bw);
% bw_positive=false(m*n,1);
% bw_positive(I_positive(positiveMapROI))=1;
% sum(bw_positive(:))
% figure; imshow(reshape(bw_positive,m,n))
   end
 end
%%
else
        %% build flag stack for succesive number NO
        flagStack_positive=false(mm,nn,NO);
        flagStack_negative=false(mm,nn,NO);
        for ii=1:NO
           file=fileName(ii).name;
           img=differentTypeRead(file,fileType);
           deltaI=img-meanImg;
            if method==2
                deltaI_I=deltaI./meanImg;
                 flagStack_positive(:,:,ii)=deltaI_I>x2;
                flagStack_negative(:,:,ii)=deltaI_I<-x2;
            elseif method==1
                 flagStack_positive(:,:,ii)=deltaI>x1*stdImg;
                flagStack_negative(:,:,ii)=deltaI<-x1*stdImg;
            elseif method==3
                flagStack_positive(:,:,ii)=img(:,:,1)>x3a;
                flagStack_negative(:,:,ii)=img(:,:,1)<x3b;
            end
        end
        hh=waitbar(0,'wait...');
        for ii=NO:p
             waitbar(ii/p,hh, ...
             [num2str(100*ii/p,'%03.1f'),'% completed']);
            file=fileName(ii).name;
            img=differentTypeRead(file,fileType);
            deltaI=img-meanImg;
               if method==2
                deltaI_I=deltaI./meanImg;
                positive1=deltaI_I>x2;
                negative1=deltaI_I<-x2;
            elseif method==1
                positive1=deltaI>x1*stdImg;
                negative1=deltaI<-x1*stdImg;
               elseif method==3
                  positive1=img(:,:,1)>x3a;
                  negative1=img(:,:,1)<x3b;
               end
               ii_2=mod(ii,NO);
               if ii_2==0
                   ii_2=NO;
               end
               flagStack_positive(:,:,ii_2)=positive1;
               flagStack_negative(:,:,ii_2)=negative1;
               positive2=flagStack_positive(:,:,1);
               for ss=2:NO
                   positive2=positive2&flagStack_positive(:,:,ss);
               end
                negative2=flagStack_negative(:,:,1);
               for ss=2:NO
                   negative2=negative2&flagStack_negative(:,:,ss);
               end
               if get(handles.cb_saveMap,'value')
                   %% save positive2 & negative2
    %                    cd(resultPath)
    %             saveName=[file(1:end-fileTypeLength),'mat'];
    %             save(saveName,'imgIOS');
    %
                    if ii==NO
                        result1='positiveMap';
                        resultPath1=fullfile(filePath,result1);
                        if exist(resultPath1)==7
                        else
                            mkdir(result1)
                        end
                        result2='negativeMap';
                        resultPath2=fullfile(filePath,result2);
                        if exist(resultPath2)==7
                        else
                            mkdir(result2)
                        end
                        result3='pnMAP';
                        resultPath3=fullfile(filePath,result3);
                        if exist(resultPath3)==7
                        else
                            mkdir(result3)
                        end
                    end
                    cd(resultPath1)
                    imwrite(positive2,[file(1:end-length(fileType)),'tif'],'tif');
                    cd(resultPath2)
                    imwrite(negative2,[file(1:end-length(fileType)),'tif'],'tif');
                    cd(resultPath3)
                    imwrite(or(positive2,negative2),[file(1:end-length(fileType)),'tif'],'tif')
                    cd(filePath)
               end
                %%
%                positive2=logical(prod(flagStack_positive,3));
%                negative2=logical(true(flagStack_negative,3));
               tk_ii=0;
                   for jj=1:length(pointPst)
                        for kk=1:length(pointPst{jj})
                            cx=pointPst{jj}{kk}(:,1);
                            cy=pointPst{jj}{kk}(:,2);
                            bw=roipoly(img,cx,cy);
                          tk_ii=tk_ii+1;
%                           positiveNO(ii-NO,tk_ii)=sum(positive2(bw));
%                           negativeNO(ii-NO,tk_ii)=sum(negative2(bw));
                          positiveNO(ii-NO+1,tk_ii)=sum(positive2(bw));
                          negativeNO(ii-NO+1,tk_ii)=sum(negative2(bw));
                          positiveSum(ii-NO+1,tk_ii)=sum(positive2(bw).*img(bw));
                          negativeSum(ii-NO+1,tk_ii)=sum(negative2(bw).*img(bw));
                        end
                   end
%                     if ii>=eval(get(handles.et_positiveNeg_from,'string')) && ...
%                         ii<=eval(get(handles.et_positiveNeg_to,'string'))
                     positiveMap=or(positiveMap,positive2);
                      negativeMap=or(negativeMap,negative2);
%                     end
        end
        positiveNO(end-(NO-1):end,:)= ones(NO,1)*positiveNO(end-(NO)+1,:);
        negativeNO(end-(NO-1):end,:)= ones(NO,1)*negativeNO(end-(NO)+1,:);
         positiveSum(end-(NO-1):end,:)= ones(NO,1)*positiveSum(end-(NO)+1,:);
        negativeSum(end-(NO-1):end,:)= ones(NO,1)*negativeSum(end-(NO)+1,:);
        close(hh)
%     elseif strcmp(fileType,'mat')~=1 %% save data as uint16
%         linesOfGroup=floor(m/zones);
%         for ii=1:zones
%             ii_2=(ii-1)*linesOfGroup+1:ii*linesOfGroup;
%             rawData=zeros(length(ii_2),n,p,'uint16');
%             for jj=1:p
%                 file=fileName(jj).name;
%                 img=imread(file);
%                 rawData(:,:,jj)=img(ii_2,:,1);
%             end
%             for kk=1:size(rawData,1)
%                 lineData=squeeze(rawData(kk,:,:));
%                 lineData2=IOS_time_gui_filter(single(lineData),handles,'onYdirection');
%
%             end
%         end
%         if mod(m,zones)~=0
%             ii_2=zones*linesOfGroup+1:m;
%             for jj=1:p
%                 file=fileName(jj).name;
%                 img=imread(file);
%                 rawData(:,:,jj)=img(ii_2,:,1);
%             end
%             for kk=1:size(rawData,1)
%                 lineData=squeeze(rawData(kk,:,:));
%                 lineData2=IOS_time_gui_filter(single(lineData),handles,'onYdirection');
%
%             end
%         end
end
 pixel_number=zeros(1,sum(ROI_N(:)));
 tt=1;
 handles.positiveMap=logical(positiveMap);
 handles.negativeMap=logical(negativeMap);
 for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        bw=roipoly(img,cx,cy);
        numberTmp=sum(bw(:));
        pixel_number(tt)=numberTmp;
        positiveNO(:,tt)=positiveNO(:,tt)/numberTmp;
        negativeNO(:,tt)=negativeNO(:,tt)/numberTmp;
        tt=tt+1;
    end
 end
 handles.pixel_number=pixel_number;
handles.positiveNO=positiveNO;
handles.negativeNO=negativeNO;
handles. positiveNegativeNO=positiveNO+negativeNO;
handles.positiveSum=positiveSum;
handles.negativeSum=negativeSum;
handles. positiveNegativeSum=positiveSum-negativeSum;
guidata(hObject,handles)

%%
% --- Executes on selection change in pm_method.
function pm_method_Callback(hObject, eventdata, handles)
% hObject    handle to pm_method (see GCBO)
handles=IOS_time_parameterRead(handles);
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns pm_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_method
function pm_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_method (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_positiveNegativeFolder_Callback(hObject, eventdata, handles)
% hObject    handle to et_positiveNegativeFolder (see GCBO)
% et_positiveNegativeFolder as text
%        str2double(get(hObject,'String')) returns contents of et_positiveNegativeFolder as a double
function et_positiveNegativeFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_positiveNegativeFolder (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton33_Callback(hObject, eventdata, handles)
positiveNO=handles.positiveNO;
negativeNO=handles.negativeNO;
positiveSum=handles.positiveSum;
negativeSum=handles.negativeSum;
positiveNegativeSum=handles. positiveNegativeSum;
positiveNegativeNO=handles. positiveNegativeNO;
p=size(positiveNO,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
colorTmp=[1,0,0;
          0,1,0;
          0,0,1;
          0.502,0,1;
          0,1,1;
          1,0,1;
          0.3,0.3,0.7;
          0.5,1/2,1/2];
pointPst=handles.pointPst;
%%
figure(50);
h1=subplot(1,3,1);
h2=subplot(1,3,2);
h3=subplot(1,3,3);
set(h1,'nextPlot','replace')
set(h2,'nextPlot','replace')
set(h3,'nextPlot','replace')
legendName=cell(size(positiveNO,2)-1,1);
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=0;
pixel_number=handles.pixel_number;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        plot(h1,timeCourse,positiveNO(:,vv),ROIii{jj},'color',colorTmp(ii,:))
         plot(h2,timeCourse,positiveSum(:,vv),ROIii{jj},'color',colorTmp(ii,:))
         ttt=positiveNO(:,vv);ttt=ttt*pixel_number(vv);
         ttt(ttt==0)=1;
         plot(h3,timeCourse,positiveSum(:,vv)./ttt,ROIii{jj},'color',colorTmp(ii,:))
        if vv==1
            set(h1,'nextPlot','add')
            set(h2,'nextPlot','add')
           % hold on;
        end
        legendName{vv}=['group',num2str(ii),'ROI',num2str(jj),'pixels',num2str(pixel_number(vv))];
    end
end
legend(h1,legendName)
legend(h2,legendName)
legend(h3,legendName)
grid on;
xlabel(timeCourseUnit)
hold off;
title(h1,'positive quantity')
title(h2,'positive intensity sum')
title(h3,'posiitve intensity mean')
%%
%%
figure(51);
h1=subplot(1,3,1);
h2=subplot(1,3,2);
h3=subplot(1,3,3);
set(h1,'nextPlot','replace')
set(h2,'nextPlot','replace')
set(h3,'nextPlot','replace')
legendName=cell(size(positiveNO,2)-1,1);
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=0;
pixel_number=handles.pixel_number;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        plot(h1,timeCourse,negativeNO(:,vv),ROIii{jj},'color',colorTmp(ii,:))
         plot(h2,timeCourse,negativeSum(:,vv),ROIii{jj},'color',colorTmp(ii,:))
         ttt=negativeNO(:,vv);ttt=ttt*pixel_number(vv);
         ttt(ttt==0)=1;
         plot(h3,timeCourse,negativeSum(:,vv)./ttt,ROIii{jj},'color',colorTmp(ii,:))
        if vv==1
            set(h1,'nextPlot','add')
            set(h2,'nextPlot','add')
           % hold on;
        end
        legendName{vv}=['group',num2str(ii),'ROI',num2str(jj),'pixels',num2str(pixel_number(vv))];
    end
end
legend(h1,legendName)
legend(h2,legendName)
legend(h3,legendName)
grid on;
xlabel(timeCourseUnit)
hold off;
title(h1,'negative quantity')
title(h2,'negative intensity sum')
title(h3,'negative intensity mean')
%%
figure(52);
h1=subplot(1,3,1);
h2=subplot(1,3,2);
h3=subplot(1,3,3);
set(h1,'nextPlot','replace')
set(h2,'nextPlot','replace')
set(h3,'nextPlot','replace')
legendName=cell(size(positiveNO,2)-1,1);
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=0;
pixel_number=handles.pixel_number;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        plot(h1,timeCourse,positiveNegativeNO(:,vv),ROIii{jj},'color',colorTmp(ii,:))
         plot(h2,timeCourse,positiveNegativeSum(:,vv),ROIii{jj},'color',colorTmp(ii,:))
         ttt=positiveNegativeNO(:,vv);ttt=ttt*pixel_number(vv);
         ttt(ttt==0)=1;
         plot(h3,timeCourse,positiveNegativeSum(:,vv)./ttt,ROIii{jj},'color',colorTmp(ii,:))
        if vv==1
            set(h1,'nextPlot','add')
            set(h2,'nextPlot','add')
           % hold on;
        end
        legendName{vv}=['group',num2str(ii),'ROI',num2str(jj),'pixels',num2str(pixel_number(vv))];
    end
end
legend(h1,legendName)
legend(h2,legendName)
legend(h3,legendName)
grid on;
xlabel(timeCourseUnit)
hold off;
title(h1,'positiveNegative quantity')
title(h2,'positiveNegative intensity sum')
title(h3,'positiveNegative intensity mean')
guidata(hObject,handles)
function pushbutton34_Callback(hObject, eventdata, handles)
cd(handles.filePath)
result=get(handles.et_positiveNegativeFolder,'string');
resultFolder=fullfile(handles.filePath,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
[file,filePath]=uiputfile('positive_negative.mat');
cd(filePath)
timeCourse=handles.timeCourse;
positiveNO=handles.positiveNO;
negativeNO=handles.negativeNO;
positiveNegativeNO=handles.positiveNegativeNO;
positiveSum=handles.positiveSum;
negativeSum=handles.negativeSum;
positiveNegativeSum=handles. positiveNegativeSum;
positiveSum=double(positiveSum);
negativeSum=double(negativeSum);
positiveNegativeSum=double(positiveNegativeSum);
pixel_number=double(handles.pixel_number);
save('pixel_number.txt','pixel_number','-ASCII');
save('positiveSum.txt','positiveSum','-ASCII');
save('negativeSum.txt','negativeSum','-ASCII');
save('positiveNegativeSum.txt','positiveNegativeSum','-ASCII');
[m,n]=size(positiveNO);
positiveSave=zeros(m,n+1);
positiveSave(:,1)=timeCourse;
positiveSave(:,2:end)=positiveNO;
save('positiveQuantity.txt','positiveSave','-ASCII');
p=length(pixel_number);
positiveMean=zeros(m,n+1);
positiveMean(:,1)=timeCourse;
for ii=1:p
    ttt=positiveNO(:,ii);ttt=ttt*pixel_number(ii);
    ttt(ttt==0)=1;
    ttt=double(ttt);
    positiveMean(:,ii+1)=positiveSum(:,ii)./ttt;
end
save('positiveMean.txt','positiveMean','-ASCII');
negativeMean=zeros(m,n+1);
negativeMean(:,1)=timeCourse;
for ii=1:p
    ttt=negativeNO(:,ii);ttt=ttt*pixel_number(ii);
    ttt(ttt==0)=1;
    ttt=double(ttt);
    negativeMean(:,ii+1)=negativeSum(:,ii)./ttt;
end
save('negativeMean.txt','negativeMean','-ASCII');
positiveNegativeMean=zeros(m,n+1);
positiveNegativeMean(:,1)=timeCourse;
for ii=1:p
    ttt=positiveNegativeNO(:,ii);ttt=ttt*pixel_number(ii);
    ttt(ttt==0)=1;
    ttt=double(ttt);
    positiveNegativeMean(:,ii+1)=positiveNegativeSum(:,ii)./ttt;
end
save('positiveNegativeMean.txt','positiveNegativeMean','-ASCII');
negativeSave=zeros(m,n+1);
negativeSave(:,1)=timeCourse;
negativeSave(:,2:end)=negativeNO;
save('negativeQuantity.txt','negativeSave','-ASCII');
positiveNegativeSave=zeros(m,n+1);
positiveNegativeSave(:,1)=timeCourse;
positiveNegativeSave(:,2:end)=positiveNegativeNO;
save('positiveNegativeQuantity.txt','positiveNegativeSave','-ASCII');
pointPst=handles.pointPst;
save(file,'positiveNO','negativeNO','positiveNegativeNO','pointPst', ...
    'positiveSum','negativeSum','positiveNegativeSum','pixel_number')
figure(50)
saveas(gcf,'positive','tiff')
figure(51)
saveas(gcf,'negative','tiff')
figure(52)
saveas(gcf,'positiveMinusNegative','tiff')
guidata(hObject,handles)
function pushbutton35_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.mat');
cd(filePath)
load(file)
handles.positiveNO=positiveNO;
handles.negativeNO=negativeNO;
handles.positiveNegativeNO=positiveNegativeNO;
handles.pointPst=pointPst;
handles.positiveSum=positiveSum;
handles.negativeSum=negativeSum;
handles. positiveNegativeSum=positiveNegativeSum;
handles.pixel_number=pixel_number;
guidata(hObject,handles)
function et_x2_Callback(hObject, eventdata, handles)
% hObject    handle to et_x2 (see GCBO)
% et_x2 as text
%        str2double(get(hObject,'String')) returns contents of et_x2 as a double
function et_x2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_x2 (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_x1_Callback(hObject, eventdata, handles)
% hObject    handle to et_x1 (see GCBO)
% et_x1 as text
%        str2double(get(hObject,'String')) returns contents of et_x1 as a double
function et_x1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_x1 (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on mouse press over axes background.
function axes_rawImg_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_rawImg (see GCBO)
pt=get(gca,'currentpoint');
pt=round(pt);
yPos=(pt(1,2)); xPos=(pt(1,1));
handles.pointer=[xPos,yPos];
handles=IOS_time_plotCurve(hObject,handles);
guidata(gca,handles)
function handles=IOS_time_plotCurve(hObject,handles)
axes(handles.axes_rawImg)
% set(gcf,'currentAxes',handles.axes_rawImg)
if isfield(handles,'crossLine')
    if ishandle(handles.crossLine)
        delete(handles.crossLine)
    end
end
hold on
a=get(gca,'xlim');
b=get(gca,'ylim');
h1=line(a,[handles.pointer(2),handles.pointer(2)]);
h2=line([handles.pointer(1),handles.pointer(1)],b);
set(h1,'color',[1,0,0])
set(h2,'color',[1,0,0])
set(allchild(gca),'ButtonDownFcn','IOS_Software(''axes_rawImg_ButtonDownFcn'',gco,[],guidata(gcbo))')
set(allchild(gca),'hitTest','on')
%%
set(handles.st_position,'string',['x=',num2str(handles.pointer(1)),'; y=', ...
    num2str(handles.pointer(2))])
t=handles.img(handles.pointer(2),handles.pointer(1),1);
t2=handles.imgRaw(handles.pointer(2),handles.pointer(1),1);
        set(handles.st_position,'string',['x=',num2str(handles.pointer(1)),'; y=', ...
    num2str(handles.pointer(2)),'; new=',num2str(t),';old=',num2str(t2)])
%%
%     h_img=findobj(handles.axes_rawImg,'type','image');
% %     h_img=differentTypeRead(fullfile(handles.filePath,handles.file),handles.fileType);
%     %set(h_img,'CData',uint16(img))
%     if isempty(h_img)==0
%         img=get(h_img,'CData');
%         t=img(handles.pointer(2),handles.pointer(1));
%         if max(img(:))<1.5
%             t=t*255;
%         end
% %         cd(handles.filePath)
%         img=differentTypeRead(fullfile(handles.filePath,handles.file),handles.fileType);
%         t2=img(handles.pointer(2),handles.pointer(1),1);
%         set(handles.st_position,'string',['x=',num2str(handles.pointer(1)),'; y=', ...
%     num2str(handles.pointer(2)),'; new=',num2str(t),';old=',num2str(t2)])
%     end
axes(handles.axes_rawImg);
% set(gcf,'currentAxes',handles.axes_rawImg)
handles.crossLine=[h1,h2];
% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%   Key: name of the key that was pressed, in lower case
%   Character: character interpretation of the key(s) that was pressed
%   Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
curr_char = int8(get(gcf,'CurrentCharacter'));
if isempty(curr_char)
    return;
end
xpos = handles.pointer(1);
ypos = handles.pointer(2);
% Keys:
% - up:   30
% - down:   31
% - left:   28
% - right:   29
% - '1': 49
% - '2': 50
% - '3': 51
% - 'e': 101
% - plus:  43
% - minus:  45
    yStep=1;
    xStep=1;
    switch curr_char
        case 30
            ypos=ypos-yStep;
        case 31
            ypos=ypos+yStep;
        case 28
            xpos=xpos-xStep;
        case 29
            xpos=xpos+xStep;
    end
handles.pointer=[xpos,ypos];
handles=IOS_time_plotCurve(hObject,handles);
guidata(hObject, handles);
% --- Executes on button press in cb_ROIsmall.
function cb_ROIsmall_Callback(hObject, eventdata, handles)
% hObject    handle to cb_ROIsmall (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of cb_ROIsmall
function pushbutton36_Callback(hObject, eventdata, handles)
function pushbutton37_Callback(hObject, eventdata, handles)
function pushbutton38_Callback(hObject, eventdata, handles)
function pushbutton39_Callback(hObject, eventdata, handles)
function et_maxNumber_Callback(hObject, eventdata, handles)
% hObject    handle to et_maxNumber (see GCBO)
% et_maxNumber as text
%        str2double(get(hObject,'String')) returns contents of et_maxNumber as a double
function et_maxNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_maxNumber (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton40_Callback(hObject, eventdata, handles)
function pushbutton41_Callback(hObject, eventdata, handles)
function pushbutton42_Callback(hObject, eventdata, handles)
function pushbutton43_Callback(hObject, eventdata, handles)
function pushbutton44_Callback(hObject, eventdata, handles)
function et_fileName_path_Callback(hObject, eventdata, handles)
% hObject    handle to et_fileName_path (see GCBO)
% et_fileName_path as text
%        str2double(get(hObject,'String')) returns contents of et_fileName_path as a double
function et_fileName_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_fileName_path (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_fileName_file_Callback(hObject, eventdata, handles)
% hObject    handle to et_fileName_file (see GCBO)
% et_fileName_file as text
%        str2double(get(hObject,'String')) returns contents of et_fileName_file as a double
function et_fileName_file_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_fileName_start_Callback(hObject, eventdata, handles)
function et_fileName_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_fileName_to_Callback(hObject, eventdata, handles)
function et_fileName_to_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_fileName_format_Callback(hObject, eventdata, handles)
function et_fileName_format_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton45_Callback(hObject, eventdata, handles)
method='userInput';
handles=imgLoadInitial(hObject,handles,method);
guidata(hObject,handles);
% --- Executes on button press in cb_rotation.

function pushbutton52_Callback(hObject, eventdata, handles)
if get(handles.et_meanFirstROI,'value')~=1
    [file,filePath]=uiputfile('referenceSignal.txt');
    cd(filePath)
    referenceSignal=double(handles.referenceSignal);
    save(file,'referenceSignal','-ASCII');
end
function pushbutton53_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
rawData=zeros(p,length(pointPst)+1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
h_wait=waitbar(0,'wait');
groupQuantity=zeros(length(pointPst),1);
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    img=differentTypeRead(file,fileType);
 if get(handles.pp_spfilterMethod,'value')>1
    img=differentTypeReadFilter_handles(img,handles);
 end    
   % vv=1;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
           % vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            ROItmp=img(bw);
            rawData(kk,ii+1)=rawData(kk,ii+1)+sum(ROItmp(:));
            if kk==1
                groupQuantity(ii)=groupQuantity(ii)+sum(bw(:));
            end
        end
        rawData(kk,ii+1)=rawData(kk,ii+1)/groupQuantity(ii);
    end
end
close(h_wait)
rawData(:,2:end)= ...
    IOS_time_gui_filter(rawData(:,2:end),handles,'vertical');
handles.rawData_meanOfgroup=rawData;
guidata(hObject,handles)
function pushbutton54_Callback(hObject, eventdata, handles)
cd(handles.filePath)
result=['meanGroup'];
% if get(handles.pp_spfilterMethod,'value')>1
%     result=[result,handles.filterMethod.saveName];
% end
resultFolder=fullfile(handles.filePath,result);
% if exist(resultFolder)==7
% else
    mkdir(result)
% end
cd(resultFolder)
[file,filePath]=uiputfile('meanOfGroup_total.txt');
cd(filePath)
figure(105)
saveas(gcf,'meanOfGroup_total','tiff')
% save(file,'meanOfGroup_total','-ASCII');
meanOfROI=handles.rawData_meanOfgroup;
rawData_meanOfgroup_IOS=handles.rawData_meanOfgroup_IOS;
meanOfROI(:,1)=rawData_meanOfgroup_IOS(:,1);
save(file,'meanOfROI','-ASCII');
save(['IOS_',file],'rawData_meanOfgroup_IOS','-ASCII');
function pushbutton55_Callback(hObject, eventdata, handles)
rawData=handles.rawData_meanOfgroup;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(104);
title('mean of group')
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-*','--*',':','-.','*-','*--','*:','*-.'};
vv=1;
for ii=1:length(pointPst)
    for jj=1:1
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'LineWidth',1,'MarkerSize',3,'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
       ccc=(rawData(:,vv));snrCCC=mean(ccc)/std(ccc);
        legendName{vv-1}=['group',num2str(ii),'SNR:',num2str(snrCCC)];
    end
end
legend(legendName)
h_legend=legend(legendName);
set(h_legend,'FontSize',16);
set(gca,'fontsize',15)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;
handles.rawData_meanOfgroup=rawData;
xlabel(timeCourseUnit)
hold off;
guidata(hObject,handles)
pre=eval(get(handles.et_pre2,'string'));
rawData(:,2:end)=dIOSization(rawData(:,2:end),pre);
%%

figure(105);
title('mean of group IOS')
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-*','--*',':','-.','*-','*--','*:','*-.'};
vv=1;
step=0;
timeToHalf=zeros(length(pointPst),1);
maxValue=zeros(length(pointPst),1);
rawData_fit=rawData;
for ii=1:length(pointPst)
    for jj=1:1
        vv=vv+1;
        if vv==2
            hold on;
        end
%     [timeToHalf(vv-1),maxValue(vv-1),rawData_fit(:,vv)]= ...
%         findHalf_expo_MLA(timeCourse,rawData(:,vv),pre,p);           
        plot(rawData(:,1),rawData(:,vv)+step*vv,ROIii{jj},'LineWidth',1,'MarkerSize',3,'color',colorTmp(ii,:))
        legendName{vv-1}=['group',num2str(ii)];


    
%        %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
% %        ccc=(rawData(:,vv));snrCCC=mean(ccc)/std(ccc);
% %         legendName{vv-1}=['group',num2str(ii),'halfT:',num2str(timeToHalf(vv-1)), ...
% %             'Max:',num2str(maxValue(vv-1))];
% %                     plot(rawData(:,1),rawData_fit(:,vv)+step*vv)
    end
end
legend(legendName)
h_legend=legend(legendName);
set(h_legend,'FontSize',16);
set(gca,'fontsize',15)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;
% xlabel(timeCourseUnit)
% hold off;
handles.rawData_meanOfgroup_IOS=rawData;
% handles.rawData_meanOfgroup_IOS_fit=rawData_fit;
guidata(hObject,handles)

function data=dIOSization(data,pre)
[m,n]=size(data);

for ii=1:n
    a=data(:,ii);
    aMean=mean(a(1:pre));
    if aMean~=0
        a=(a-aMean)/aMean;
        data(:,ii)=a;
    
    end
end

function data=dIization(data,pre)
[m,n]=size(data);

for ii=1:n
    a=data(:,ii);
    aMean=mean(a(1:pre));
    if aMean~=0
          a=(a-aMean);
        data(:,ii)=a;
    
    end
end

function [timeToHalf,maxValue,y4]=findHalf_expo_MLA(x0,y0,n1,n2,varargin)
baseY=mean(y0(1:n1));
if abs(x0(n1))>0.0001
    error('time scale is not right')
end
x=x0(n1:n2);
x=x(:);
y=y0(n1:n2);
y=y(:);
fun=inline('y-baseY+par(1)-par(1)*exp(-x/par(2))','par','y','x','baseY');
fun2=inline('baseY-par(1)+par(1)*exp(-x/par(2))','par','x','baseY');
par0=[0.01    .1];
[par,resnorm,residual]=lsqnonlin(fun,par0,[],[],[],y,x,baseY);
timeToHalf=-par(2)*log(0.5);
% figure(3);hold off; plot(x,fun2(par,x,baseY),x,y,'*')

    y2=fun2(par,x,baseY);
    y3=fun2(par,zeros(n1-1,1),baseY);
    y4=y0(1:n2);y4(1:n1-1)=y3;y4(n1:end)=y2;
    
%     figure(4);hold off; plot(x0,y0,'+');grid on;
%     hold on; plot(x0(1:n2),y4,'r-','lineWidth',2)
%     title(['halfTime:',num2str(timeToHalf)])

maxValue=baseY-par(1);
     

function pushbutton56_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.txt');
cd(filePath)
handles.rawData_meanOfgroup=load(file);
guidata(hObject,handles)
% --- Executes on button press in pb_tab_IOS.
function pb_tab_IOS_Callback(hObject, eventdata, handles)
% hObject    handle to pb_tab_IOS (see GCBO)
% --- Executes on button press in pb_DIOS.
function pb_DIOS_Callback(hObject, eventdata, handles)
% hObject    handle to pb_DIOS (see GCBO)
function pushbutton60_Callback(hObject, eventdata, handles)
% --- Executes on button press in pb_DIOS_1.
function pb_DIOS_1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_DIOS_1 (see GCBO)
% --- Executes on button press in pb_vertical_tab.
function pb_vertical_tab_Callback(hObject, eventdata, handles)
% hObject    handle to pb_vertical_tab (see GCBO)
handles.verticalTab=[handles.imgQuality,handles.cb_ROIuserInput, ...
    handles.t_positions,handles.pushbutton30,handles.et_verticalLine, ...
    handles.text28, handles.et_fig,handles.pb_verticalTab1];
handles.DIOStab=[handles.pb_DIOS_1,handles.st_dIOS_to,handles.et_dIOS_to,...
    handles.st_dIOS_folderName, handles.et_dIOS_folderName, handles.et_dIOS_conversion, ...
    handles.pb_dIOS_mat, handles.pb_dIOS_single, handles.pb_dIOS_color, ...
    handles.st_dIOS_scale,handles.st_dIOS_binning, handles.et_dIOS_binning, ...
    handles.pb_dIOS_cover2, handles.et_dIOS_pseudocolor, handles.et_dIOS_gray, ...
    handles.text62,handles.pp_dIOScolormap,handles.text60,handles.cb_inverseNegative_dIOS,handles.cb_continuous];
set(handles.verticalTab,'visible','on')
set(handles.DIOStab,'visible','off')
set(handles.pb_dIOS_tab,'foregroundColor',[1,0,0])
set(handles.pb_vertical_tab,'foregroundColor',[0,0,0])
guidata(hObject,handles)
% --- Executes on button press in pb_dIOS_tab.
function pb_dIOS_tab_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dIOS_tab (see GCBO)
handles.verticalTab=[handles.imgQuality,handles.cb_ROIuserInput, ...
    handles.t_positions,handles.pushbutton30,handles.et_verticalLine, ...
    handles.text28, handles.et_fig,handles.pb_verticalTab1];
handles.DIOStab=[handles.pb_DIOS_1,handles.st_dIOS_to,handles.et_dIOS_to,...
    handles.st_dIOS_folderName, handles.et_dIOS_folderName, handles.et_dIOS_conversion, ...
    handles.pb_dIOS_mat, handles.pb_dIOS_single, handles.pb_dIOS_color, ...
    handles.st_dIOS_scale,handles.st_dIOS_binning, handles.et_dIOS_binning, ...
    handles.pb_dIOS_cover2, handles.et_dIOS_pseudocolor, handles.et_dIOS_gray, ...
    handles.text62,handles.pp_dIOScolormap,handles.text60,handles.cb_inverseNegative_dIOS,handles.cb_continuous];
set(handles.verticalTab,'visible','off')
set(handles.DIOStab,'visible','on')
set(handles.pb_dIOS_tab,'foregroundColor',[0,0,0])
set(handles.pb_vertical_tab,'foregroundColor',[1,0,0])
guidata(hObject,handles)
% --- Executes on button press in pb_verticalTab1.
function pb_verticalTab1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_verticalTab1 (see GCBO)
function et_dIOS_to_Callback(hObject, eventdata, handles)
% hObject    handle to et_dIOS_to (see GCBO)
% et_dIOS_to as text
%        str2double(get(hObject,'String')) returns contents of et_dIOS_to as a double
function et_dIOS_to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_dIOS_to (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_dIOS_mat.
function pb_dIOS_mat_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dIOS_mat (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
baselineNumber=eval(get(handles.et_dIOS_to,'string'));
binning=eval(get(handles.et_dIOS_binning,'string'));
p=length(fileName);
fileTypeLength=length(fileType);
%% get result directory
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1=get(handles.et_dIOS_folderName,'string');
resultPath1=fullfile(filePathTmp,resultName1);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
cd(filePath)
img2=differentTypeRead(file,fileType);
[mm,nn,nn_trash]=size(img2);
inverseNegativeFlag=get(handles.cb_inverseNegative_dIOS,'value');
if get(handles.cb_continuous,'value')
        imgAvgStack=zeros(mm,nn,baselineNumber);
%         imgSingleStack=zeros(mm,nn,binning);
        %% get the first baselineNumber image's average
        for ii=1:baselineNumber
            file=fileName(ii).name;
                        img2=differentTypeRead(file,handles.fileType);
            imgAvgStack(:,:,ii)=(img2(:,:,1));
        end
        %% binning effect
        n=p;
        h_wait=waitbar(0,'wait');
        ii=baselineNumber;
        file=fileName(ii).name;
        img2=differentTypeRead(file,handles.fileType);
        imgSingleStack2=img2(:,:,1);
        for ii=baselineNumber+1:n-binning+1
                  waitbar(ii/n,h_wait, ...
             [num2str(100*ii/n,'%03.1f'),'% completed']);
            imgSingleStack=zeros(mm,nn,binning);
            cd(filePath)
             for jj=ii:ii+binning-1
                 ss=jj-ii+1;
                file=fileName(jj).name;
                img2=differentTypeRead(file,handles.fileType);
                imgSingleStack(:,:,ss)=(img2(:,:,1));
             end
                imgAvgStack(:,:,1:end-1)=imgAvgStack(:,:,1+1:end);
                imgAvgStack(:,:,end-1+1:end)=imgSingleStack2;
            file=fileName(ii).name;
            img2=differentTypeRead(file,handles.fileType);
            imgSingleStack2=img2(:,:,1);
            imgSingle=single(mean(imgSingleStack,3));
            imgAvg=single(mean(imgAvgStack,3));
            nanFlag=imgAvg==0;
            if sum(nanFlag(:))~=0
                imgAvg=imgAvg+1;
            end
            imgSave=(imgSingle-imgAvg)./imgAvg;
            if inverseNegativeFlag
                imgSave=abs(imgSave);
            end
            cd(resultPath1)
            saveName=[file(1:end-fileTypeLength),'mat'];
            save(saveName,'imgSave');
        end
        close(h_wait)
else
    if strcmp(fileType,'mat')~=1 && binning<baselineNumber
        imgAvgStack=zeros(mm,nn,baselineNumber,'uint16');
        %% get the first baselineNumber image's average
        for ii=1:baselineNumber
            file=fileName(ii).name;
            img2=imread(file);
            imgAvgStack(:,:,ii)=uint16(img2(:,:,1));
        end
        %% binning effect
        n=floor((p-baselineNumber)/binning);
        h_wait=waitbar(0,'wait');
        for ii=1:n
                  waitbar(ii/n,h_wait, ...
             [num2str(100*ii/n,'%03.1f'),'% completed']);
            imgSingleStack=zeros(mm,nn,binning,'uint16');
            cd(filePath)
             for jj=baselineNumber+binning*(ii-1)+1:baselineNumber+binning*(ii-1)+binning
                 ss=jj-(baselineNumber+binning*(ii-1));
                file=fileName(jj).name;
                img2=imread(file);
                imgSingleStack(:,:,ss)=uint16(img2(:,:,1));
             end
            if ii~=1
                imgAvgStack(:,:,1:end-binning)=imgAvgStack(:,:,1+binning:end);
                imgAvgStack(:,:,end-binning+1:end)=imgSingleStack2;
            end
            imgSingleStack2=imgSingleStack;
            imgSingle=single(mean(imgSingleStack,3));
            imgAvg=single(mean(imgAvgStack,3));
            nanFlag=imgAvg==0;
            if sum(nanFlag(:))~=0
                imgAvg=imgAvg+1;
            end
            imgSave=(imgSingle-imgAvg)./imgAvg;
            if inverseNegativeFlag
                imgSave=abs(imgSave);
            end
            cd(resultPath1)
            saveName=[fileName(baselineNumber+binning*(ii-1)+1).name(1:end-fileTypeLength),'mat'];
            save(saveName,'imgSave');
        end
        close(h_wait)
    else
        imgAvgStack=zeros(mm,nn,baselineNumber);
        %% get the first baselineNumber image's average
        for ii=1:baselineNumber
            file=fileName(ii).name;
            img2=differentTypeRead(file,handles.fileType);
            imgAvgStack(:,:,ii)=(img2(:,:,1));
        end
        %% binning effect
        n=floor((p-baselineNumber)/binning);
        h_wait=waitbar(0,'wait');
        for ii=1:n
                  waitbar(ii/n,h_wait, ...
             [num2str(100*ii/n,'%03.1f'),'% completed']);
            imgSingleStack=zeros(mm,nn,binning);
            cd(filePath)
             for jj=baselineNumber+binning*(ii-1)+1:baselineNumber+binning*(ii-1)+binning
                 ss=jj-(baselineNumber+binning*(ii-1));
                file=fileName(jj).name;
                img2=differentTypeRead(file,handles.fileType);
                imgSingleStack(:,:,ss)=(img2(:,:,1));
             end
            if ii~=1
                imgAvgStack(:,:,1:end-binning)=imgAvgStack(:,:,1+binning:end);
                imgAvgStack(:,:,end-binning+1:end)=imgSingleStack2;
            end
            imgSingleStack2=imgSingleStack;
            imgSingle=single(mean(imgSingleStack,3));
            imgAvg=single(mean(imgAvgStack,3));
            nanFlag=imgAvg==0;
            if sum(nanFlag(:))~=0
                imgAvg=imgAvg+1;
            end
            imgSave=(imgSingle-imgAvg)./imgAvg;
            if inverseNegativeFlag
                imgSave=abs(imgSave);
            end
            cd(resultPath1)
            saveName=[fileName(baselineNumber+binning*(ii-1)+1).name(1:end-fileTypeLength),'mat'];
            save(saveName,'imgSave');
        end
        close(h_wait)
        %%
    %     uiwait(msgbox('codes for mat type are not optimized','Title','modal'));
    %     n=floor((p-baselineNumber)/binning);
    %      h_wait=waitbar(0,'wait');
    %     for ii=1:n
    %                waitbar(ii/n,h_wait, ...
    %          [num2str(100*ii/n,'%03.1f'),'% completed']);
    %         imgSingle=zeros(mm,nn,'single');
    %         imgAvg=imgSingle;
    %         cd(filePath)
    %          for jj=baselineNumber+binning*(ii-1)+1:baselineNumber+binning*(ii-1)+binning
    %
    %             file=fileName(jj).name;
    %             img2=differentTypeRead(file,fileType);
    %             imgSingle=imgSingle+single(img2(:,:,1));
    %          end
    %          imgSingle=imgSingle/binning;
    %
    %          for kk=binning*(ii-1)+1:baselineNumber+binning*(ii-1)
    %              file=fileName(kk).name;
    %              img2=differentTypeRead(file,fileType);
    %              imgAvg=imgAvg+single(img2(:,:,1));
    %          end
    %          imgAvg=imgAvg/baselineNumber;
    %          nanFlag=imgAvg==0;
    %         if sum(nanFlag(:))~=0
    %             imgAvg=imgAvg+1;
    %         end
    %
    %         imgSave=(imgSingle-imgAvg)./imgAvg;
    %         if inverseNegativeFlag
    %             imgSave=abs(imgSave);
    %         end
    %         cd(resultPath1)
    %         saveName=[fileName(baselineNumber+binning*(ii-1)+1).name(1:end-fileTypeLength),'mat'];
    %         save(saveName,'imgSave');
    %     end
    %     close(h_wait)
    end
end
guidata(hObject,handles)
function et_dIOS_folderName_Callback(hObject, eventdata, handles)
% hObject    handle to et_dIOS_folderName (see GCBO)
% et_dIOS_folderName as text
%        str2double(get(hObject,'String')) returns contents of et_dIOS_folderName as a double
function et_dIOS_folderName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_dIOS_folderName (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton31_CreateFcn(hObject, eventdata, handles)
% --- Executes on button press in pb_dIOS_single.
function pb_dIOS_single_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dIOS_single (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
baselineNumber=eval(get(handles.et_dIOS_to,'string'));
binning=eval(get(handles.et_dIOS_binning,'string'));
p=length(fileName);
fileTypeLength=length(fileType);
selectedFileValue=get(handles.lb_curFileNames,'Value');
cd(filePath)
img2=differentTypeRead(file,fileType);
[mm,nn,nn_trash]=size(img2);
imgSingle=zeros(mm,nn,'single');
imgAvg=imgSingle;
cd(filePath)
 for jj=selectedFileValue:selectedFileValue+binning-1
    file=fileName(jj).name;
    img2=differentTypeRead(file,fileType);
    imgSingle=imgSingle+single(img2(:,:,1));
 end
 imgSingle=imgSingle/binning;
 for kk=selectedFileValue-baselineNumber:selectedFileValue-1
     file=fileName(kk).name;
     img2=differentTypeRead(file,fileType);
     imgAvg=imgAvg+single(img2(:,:,1));
 end
 imgAvg=imgAvg/baselineNumber;
 nanFlag=imgAvg==0;
if sum(nanFlag(:))~=0
    imgAvg=imgAvg+1;
end
imgSave=(imgSingle-imgAvg)./imgAvg;
convertFunction=get(handles.et_dIOS_conversion,'string');
if isfield(handles,'ff_dIOS')
    if ishandle(handles.ff_dIOS)
        close(handles.ff_dIOS)
    end
end
handles.ff_dIOS=singleRetrieve_dIOS(imgSave,convertFunction);
guidata(hObject,handles)
% --- Executes on button press in pb_dIOS_color.
function pb_dIOS_color_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dIOS_color (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
baselineNumber=eval(get(handles.et_dIOS_to,'string'));
binning=eval(get(handles.et_dIOS_binning,'string'));
p=length(fileName);
fileTypeLength=length(fileType);
%% get result directory
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1=get(handles.et_dIOS_folderName,'string');
resultPath1=fullfile(filePathTmp,resultName1);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
%%
filePathTmp=resultPath1;
resultName1=get(handles.et_dIOS_gray,'string');
resultName2=get(handles.et_dIOS_pseudocolor,'string');
resultPath1=fullfile(filePathTmp,resultName1);
resultPath2=fullfile(filePathTmp,resultName2);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
if exist(resultPath2)==7
else
    mkdir(resultName2);
end
%% define conversion function
newIntensity=inline(get(handles.et_dIOS_conversion,'string'));
%% define colormap
contents = cellstr(get(handles.pp_dIOScolormap,'String'));
colorSelectedTmp=contents(get(handles.pp_dIOScolormap,'value'));
colorSelected=eval(colorSelectedTmp{1});
%%
cd(filePath)
img2=differentTypeRead(file,fileType);
[mm,nn,nn_trash]=size(img2);
if strcmp(fileType,'mat')~=1 && binning<baselineNumber
    imgAvgStack=zeros(mm,nn,baselineNumber,'uint16');
    %% get the first baselineNumber image's average
    for ii=1:baselineNumber
        file=fileName(ii).name;
        img2=imread(file);
        imgAvgStack(:,:,ii)=uint16(img2(:,:,1));
    end
    %% binning effect
    n=floor((p-baselineNumber)/binning);
    h_wait=waitbar(0,'wait');
    for ii=1:n
              waitbar(ii/n,h_wait, ...
         [num2str(100*ii/n,'%03.1f'),'% completed']);
        imgSingleStack=zeros(mm,nn,binning,'uint16');
        cd(filePath)
         for jj=baselineNumber+binning*(ii-1)+1:baselineNumber+binning*(ii-1)+binning
             ss=jj-(baselineNumber+binning*(ii-1));
            file=fileName(jj).name;
            img2=imread(file);
            imgSingleStack(:,:,ss)=uint16(img2(:,:,1));
         end
        if ii~=1
            imgAvgStack(:,:,1:end-binning)=imgAvgStack(:,:,1+binning:end);
            imgAvgStack(:,:,end-binning+1:end)=imgSignleStack2;
        end
        imgSignleStack2=imgSingleStack;
        imgSingle=single(mean(imgSingleStack,3));
        imgAvg=single(mean(imgAvgStack,3));
        nanFlag=imgAvg==0;
        if sum(nanFlag(:))~=0
            imgAvg=imgAvg+1;
        end
        imgSave=(imgSingle-imgAvg)./imgAvg;
        img3=uint8(newIntensity(imgSave));
        cd(resultPath1)
        saveName1=[file(1:end-length(fileType)),'tif'];
        imwrite(img3,saveName1,'tiff','compression','none');
        cd(resultPath2)
        saveName2=[file(1:end-length(fileType)),'png'];
        imgColor=ind2rgb(img3,colorSelected);
        imwrite(imgColor,saveName2,'png')
    end
    close(h_wait)
else
    uiwait(msgbox('codes for mat type are not optimized','Title','modal'));
    n=floor((p-baselineNumber)/binning);
     h_wait=waitbar(0,'wait');
    for ii=1:n
               waitbar(ii/n,h_wait, ...
         [num2str(100*ii/n,'%03.1f'),'% completed']);
        imgSingle=zeros(mm,nn,'single');
        imgAvg=imgSingle;
        cd(filePath)
         for jj=baselineNumber+binning*(ii-1)+1:baselineNumber+binning*(ii-1)+binning
            file=fileName(jj).name;
            img2=differentTypeRead(file,fileType);
            imgSingle=imgSingle+single(img2(:,:,1));
         end
         imgSingle=imgSingle/binning;
         for kk=binning*(ii-1)+1:baselineNumber+binning*(ii-1)
             file=fileName(kk).name;
             img2=differentTypeRead(file,fileType);
             imgAvg=imgAvg+single(img2(:,:,1));
         end
         imgAvg=imgAvg/baselineNumber;
         nanFlag=imgAvg==0;
        if sum(nanFlag(:))~=0
            imgAvg=imgAvg+1;
        end
        imgSave=(imgSingle-imgAvg)./imgAvg;
         img3=uint8(newIntensity(imgSave));
        cd(resultPath1)
        saveName1=[file(1:end-length(fileType)),'tif'];
        imwrite(img3,saveName1,'tiff','compression','none');
        cd(resultPath2)
        saveName2=[file(1:end-length(fileType)),'png'];
        imgColor=ind2rgb(img3,colorSelected);
        imwrite(imgColor,saveName2,'png')
    end
    close(h_wait)
end
guidata(hObject,handles)
function et_dIOS_conversion_Callback(hObject, eventdata, handles)
% hObject    handle to et_dIOS_conversion (see GCBO)
% et_dIOS_conversion as text
%        str2double(get(hObject,'String')) returns contents of et_dIOS_conversion as a double
function et_dIOS_conversion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_dIOS_conversion (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_dIOS_binning_Callback(hObject, eventdata, handles)
function et_dIOS_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pm_colormap_raw_Callback(hObject, eventdata, handles)
lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles)
function pm_colormap_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_x3a_Callback(hObject, eventdata, handles)
function et_x3a_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_x3b_Callback(hObject, eventdata, handles)
function et_x3b_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_dIOS_gray_Callback(hObject, eventdata, handles)
function et_dIOS_gray_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_dIOS_pseudocolor_Callback(hObject, eventdata, handles)
function et_dIOS_pseudocolor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pb_dIOS_cover2_Callback(hObject, eventdata, handles)
function pp_dIOScolormap_Callback(hObject, eventdata, handles)
function pp_dIOScolormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_positiveNeg_from_Callback(hObject, eventdata, handles)
function et_positiveNeg_from_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function et_positiveNeg_to_Callback(hObject, eventdata, handles)
function et_positiveNeg_to_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton78_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.mat');
load(fullfile(filePath,file))
handles.rawData_positive_meanOfROI=rawData_positive_meanOfROI;
handles.rawData_pn_meanOfROI=rawData_pn_meanOfROI;
handles.rawData_negative_meanOfROI=rawData_negative_meanOfROI;
handles.pointPst=pointPst;
handles.positiveMap=positiveMap;
handles.negativeMap=negativeMap;
handles.colorTmp=colorTmp;
guidata(hObject,handles)
function handles=pushbutton79_Callback(hObject, eventdata, handles)
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
pnYellow=and(positiveMap,handles.negativeMap);
positiveMap2=xor(pnYellow,positiveMap);
negativeMap2=xor(pnYellow,negativeMap);
%% positive
rawData=handles.rawData_positive_meanOfROI;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;

legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
%% draw all together
if size(rawData,2)==2
    pre=eval(get(handles.et_pre2,'string'));
    meanP=mean(rawData(1:pre,2));

    pPlot=(rawData(:,2)-meanP)/meanP;
    
    nPlot=handles.rawData_negative_meanOfROI(:,2);
    meanN=mean(nPlot(1:pre));
    nPlot=(nPlot-meanN)/meanN;
    
    
    cx=pointPst{1}{1}(:,1);
    cy=pointPst{1}{1}(:,2);
    bwTmp=roipoly(double(positiveMap),cx,cy);
    bwTmpNumber=and(bwTmp,positiveMap2);
    pNumber=sum(bwTmpNumber(:));
    bwTmpNumber=and(bwTmp,negativeMap2);
    nNumber=sum(bwTmpNumber(:));
    pnPlot=(rawData(:,2)*pNumber+handles.rawData_negative_meanOfROI(:,2)*nNumber)/(nNumber+pNumber);
    meanPN=mean(pnPlot(1:pre));
    pnPlot=(pnPlot-meanPN)/meanPN;
    pPlot= IOS_time_gui_filter(pPlot,handles,'vertical');
    nPlot= IOS_time_gui_filter(nPlot,handles,'vertical');
    figure(108);
    plot(timeCourse,pPlot,'r',timeCourse,nPlot,'g','lineWidth',3)
    xlabel(timeCourseUnit)
    handles.rawData_pnPlot=[timeCourse(:),pPlot(:),nPlot(:)];
    
%     plot(timeCourse,pPlot,timeCourse,nPlot,timeCourse,pnPlot)
%     bwTmpNumberSum=sum(bwTmpNumber(:));
end

%% 
pPlot=handles.rawData_positive_meanOfROI(:,2:end);
nPlot=handles.rawData_negative_meanOfROI(:,2:end);
pre=eval(get(handles.et_pre2,'string'));
for ii=1:size(pPlot,2)
    meanP=mean(pPlot(1:pre,ii));
    if abs(meanP)>0.00001
        pPlot(:,ii)=(pPlot(:,ii)-meanP)/meanP;
        
    end
    meanN=mean(nPlot(1:pre,ii));
    if abs(meanN)>0.00001
        nPlot(:,ii)=(nPlot(:,ii)-meanN)/meanN;
        
    end    
end
handles.rawData_pnPlot=[timeCourse(:),pPlot,nPlot];
figure(108);close(108);figure(108)
vv=1;
minN=min(nPlot(:));
maxP=max(pPlot(:));
maxP=max([maxP,abs(minN)]);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        subplot(1,size(pPlot,2),vv)
        vv=vv+1;
%         if vv==2
%             hold on;
%         end
        plot(timeCourse,pPlot(:,vv-1),'r-')
        hold on;
        plot(timeCourse,nPlot(:,vv-1),'b-')
        ylim([-maxP,maxP])
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,positiveMap2);
        legend(['group',num2str(ii),'ROI',num2str(jj)]);
    end
end
% legend(legendName)
%%
vv=1;
figure(105);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,positiveMap2);
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum(bwTmpNumber(:)))];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on; title('positive')
handles.rawData_positive_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;
%% negative
rawDataPositive=rawData;
rawData=handles.rawData_negative_meanOfROI;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(106);
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
pre2=eval(get(handles.et_pre2,'string'));
rawDataAbs=zeros(size(rawData));
rawDataAbs(:,1)=timeCourse;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,negativeMap2);
            
            bwTmp2=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber2=and(bwTmp2,positiveMap2);      
            sum1=sum(bwTmpNumber(:));
            sum2=sum(bwTmpNumber2(:));
            rawDataAbs(:,vv)=-sum1*(rawData(:,vv)-mean(rawData(1:pre2,vv))) ...
                +sum2*(rawDataPositive(:,vv)-mean(rawDataPositive(1:pre2,vv)));
            rawDataAbs(:,vv)=rawDataAbs(:,vv)/(sum1+sum2);
            rawDataAbs(:,vv)=rawDataAbs(:,vv)+(mean(rawData(1:pre2,vv))+ ...
                mean(rawDataPositive(1:pre2,vv)))/2;
            cc2=corrcoef(rawDataPositive(:,vv),rawData(:,vv));
            cc=cc2(1,2);
            
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum(bwTmpNumber(:))),'cc',num2str(cc)];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;title('negative')
handles.rawData_negative_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;

%% absolute
rawData=rawDataAbs;
figure(110);
vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,negativeMap2);
            
            bwTmp2=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber2=and(bwTmp2,positiveMap2);      
            sum1=sum(bwTmpNumber(:));
            sum2=sum(bwTmpNumber2(:));

        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum1+sum2)];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;title('absolute')
handles.rawData_abs_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;
%% positive and negative
rawData=handles.rawData_pn_meanOfROI;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(107);
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                        cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,pnYellow);
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum(bwTmpNumber(:)))];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on; title('positive and negative mean')
handles.rawData_pn_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;
%% map
% figure(107);
% imshow(handles.positiveMap)
%
% hold on;
% title('positive map')
% colorTmp=handles.colorTmp;
% pointPst=handles.pointPst;
% for ii=1:length(pointPst)
%     for jj=1:length(pointPst{ii})
%         cx=pointPst{ii}{jj}(:,1);
%         cy=pointPst{ii}{jj}(:,2);
%
%         hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
% %             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
% %             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
%     end
% end
% hold off;
%
% figure(108);
% imshow(handles.negativeMap)
% hold on;
% title('negative map')
% colorTmp=handles.colorTmp;
% pointPst=handles.pointPst;
% for ii=1:length(pointPst)
%     for jj=1:length(pointPst{ii})
%         cx=pointPst{ii}{jj}(:,1);
%         cy=pointPst{ii}{jj}(:,2);
%
%         hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
% %             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
% %             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
%     end
% end
% hold off;
figure(109);
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
[mm,nn,nn_trash]=size(negativeMap);
pnMap=zeros(mm,nn,3,'uint8');
pnMapTmp=zeros(mm,nn,'uint8');
pnMapTmp(positiveMap)=255;
pnMap(:,:,1)=pnMapTmp;
pnMapTmp=zeros(mm,nn,'uint8');
pnMapTmp(negativeMap)=255;
pnMap(:,:,2)=pnMapTmp;
if isfield(handles,'file') && isfield(handles,'filePath')
    cd(handles.filePath)
    handles.newIntensity=inline(get(handles.et_inline,'string'));
    fileType=handles.fileType;
    file=handles.file;
    img=differentTypeRead(file,fileType);
    newImg=handles.newIntensity(img);
    newImg=uint8(newImg);
    newImg(positiveMap)=0;
    newImg(negativeMap)=0;
    pnMap(:,:,1)=pnMap(:,:,1)+newImg;
    pnMap(:,:,2)=pnMap(:,:,2)+newImg;
    pnMap(:,:,3)=pnMap(:,:,3)+newImg;
end
imshow(pnMap)
handles.pnMap=pnMap;
title('red: positive; green: negative; yellow:both')
hold on;
colorTmp=handles.colorTmp;
pointPst=handles.pointPst;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
%             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
%             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
    end
end
hold off;
guidata(hObject,handles)
function handles=pushbutton80_Callback(hObject, eventdata, handles)
cd(handles.filePath)
result=get(handles.et_positiveNegativeFolder,'string');
    if get(handles.cb_pnIOS_filter,'value')
        result=[result,handles.filterMethod2.saveName];
    end
resultFolder=fullfile(handles.filePath,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
% [file,filePath]=uiputfile('meanOfROIStatistics.mat');
% cd(filePath)
file='meanOfROIStatistics.mat';
file=file(1:end-4);
rawData_positive_meanOfROI=handles.rawData_positive_meanOfROI;
rawData_negative_meanOfROI=handles.rawData_negative_meanOfROI;
rawData_abs_meanOfROI=handles.rawData_abs_meanOfROI;
rawData_pn_meanOfROI=handles.rawData_pn_meanOfROI;
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
pointPst=handles.pointPst;
colorTmp=handles.colorTmp;
pnMap=handles.pnMap;
rawData_pnPlot=handles.rawData_pnPlot;
    positiveRawAll=handles.positiveRawAll;
    negativeRawAll= handles.negativeRawAll;
save([file,'_positive.txt'],'rawData_positive_meanOfROI','-ASCII');
save([file,'_negative.txt'],'rawData_negative_meanOfROI','-ASCII');
save([file,'_abs.txt'],'rawData_abs_meanOfROI','-ASCII');
save([file,'pn_debase.txt'],'rawData_pnPlot','-ASCII')
save([file,'_positiveNegative.txt'],'rawData_pn_meanOfROI','-ASCII');
save([file,'.mat'],'rawData_positive_meanOfROI','rawData_negative_meanOfROI', ...
    'positiveMap','negativeMap','pointPst','colorTmp', 'rawData_pn_meanOfROI','positiveRawAll','negativeRawAll')
figure(105)
saveas(gcf,'meanOfROI_positive','tiff')
figure(106)
saveas(gcf,'meanOfGroups_negative','tiff')
figure(107)
saveas(gcf,'meanOfGroups_positiveNegativeMix','tiff')
figure(109)
saveas(gcf,'distribution','tiff')
figure(108)
saveas(gcf,'PN','tiff')
figure(110)
saveas(gcf,'meanOfROI_abs','tiff')
imwrite(pnMap,'pnMap_raw.tif')
imwrite(positiveMap,'positiveMap.tif');
imwrite(negativeMap,'negativeMap.tif');
a=or(positiveMap,negativeMap);
imwrite(a,'pos_or_neg.tif')


% avg

function handles=pushbutton81_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
rawData_positive=zeros(p,sum(ROI_N(:))+1);
rawData_negative=zeros(p,sum(ROI_N(:))+1);
rawData_pn=zeros(p,sum(ROI_N(:))+1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData_positive(:,1)=timeCourse;
rawData_negative(:,1)=timeCourse;
rawData_pn(:,1)=timeCourse;
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
pnYellow=and(positiveMap,handles.negativeMap);
positiveMap2=xor(pnYellow,positiveMap);
negativeMap2=xor(pnYellow,negativeMap);
positiveRawAll=zeros(p,sum(positiveMap2(:)),'single');
negativeRawAll=zeros(p,sum(negativeMap2(:)),'single');

h_wait=waitbar(0,'wait');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    img=differentTypeRead(file,fileType);
    positiveRawAll(kk,:)=(img(positiveMap2)).';
     negativeRawAll(kk,:)=(img(negativeMap2)).';
    if get(handles.cb_pnIOS_filter,'value')
        img=differentTypeReadFilter_handles(img,handles);
    end
    
    vv=1;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
            vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(img(:,:,1),cx,cy);
            bw_positive=and(bwTmp,positiveMap2);
            bw_negative=and(bwTmp,negativeMap2);
%             bw_positive=and(bwTmp,handles.positiveMap);
%             bw_negative=and(bwTmp,handles.negativeMap);
            bw_pn=and(bwTmp,pnYellow);
            ROItmp_positive=img(bw_positive);
            ROItmp_negative=img(bw_negative);
            ROItmp_pn=img(bw_pn);
            rawData_positive(kk,vv)=mean(ROItmp_positive(:));
            rawData_negative(kk,vv)=mean(ROItmp_negative(:));
            rawData_pn(kk,vv)=mean(ROItmp_pn(:));
        end
    end
end
close(h_wait)
rawData_positive(:,2:end)= ...
    IOS_time_gui_filter(rawData_positive(:,2:end),handles,'vertical');
rawData_negative(:,2:end)= ...
    IOS_time_gui_filter(rawData_negative(:,2:end),handles,'vertical');
rawData_pn(:,2:end)= ...
    IOS_time_gui_filter(rawData_pn(:,2:end),handles,'vertical');
handles.rawData_positive_meanOfROI=rawData_positive;
handles.rawData_negative_meanOfROI=rawData_negative;
handles.rawData_pn_meanOfROI=rawData_pn;
    handles.positiveRawAll=positiveRawAll;
     handles.negativeRawAll=negativeRawAll;
guidata(hObject,handles)
function pushbutton82_Callback(hObject, eventdata, handles)
function pushbutton83_Callback(hObject, eventdata, handles)
% --- Executes on button press in cb_inverseNegative_dIOS.
function cb_inverseNegative_dIOS_Callback(hObject, eventdata, handles)
% hObject    handle to cb_inverseNegative_dIOS (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of cb_inverseNegative_dIOS
% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% --------------------------------------------------------------------
function me_imgLoad_Callback(hObject, eventdata, handles)
% hObject    handle to me_imgLoad (see GCBO)
method='imgLoad';
handles=imgLoadInitial(hObject,handles,method);
guidata(hObject,handles);
% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% --------------------------------------------------------------------
function me_texture_uniformity_Callback(hObject, eventdata, handles)
% hObject    handle to me_texture_uniformity (see GCBO)
% map function mapping the raw range to [0,1]
handles.newIntensityGroup=inline(get(handles.et_inline_texture,'string'));
% get the bins.
edgeBin=eval(get(handles.et_texture_edge,'string'));
%pBin=length(edgeBin);
centerBin=edgeBin(1:end-1)+edgeBin(2:end);
centerBin=centerBin/2;
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
rawData=cell(groupN,1);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        bw=roipoly(img(:,:,1),cx,cy);
        ROI_N(ii)=ROI_N(ii)+sum(bw(:));
    end
    rawData{ii}=zeros(ROI_N(ii),1);
end
groupQuantity=ROI_N;
uniformity=zeros(p,groupN+1);
smoothness=zeros(p,groupN+1);
% time course
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
uniformity(:,1)=timeCourse;
smoothness(:,1)=timeCourse;
h_wait=waitbar(0,'wait');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    img=differentTypeRead(file,fileType);
   % vv=1;
    for ii=1:length(pointPst)
        startIndex=1;
        for jj=1:length(pointPst{ii})
           % vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            ROItmp=img(bw);
            ROI_length=length(ROItmp);
            rawData{ii}(startIndex:startIndex+ROI_length-1)=ROItmp;
            startIndex=startIndex+ROI_length;
        end
        rawData{ii}=handles.newIntensityGroup(rawData{ii});
        % histogram
        yBinTmp=histc(rawData{ii}(:),edgeBin);
        %normalized
        yBinTmp=yBinTmp/groupQuantity(ii);
        % add the bin on the edge to the last
        yBinTmp(end-1)=yBinTmp(end-1)+yBinTmp(end);
        yBin=yBinTmp(1:end-1);
 %       percentTotal=sum(yBin(:))/mm/nn;
 %% uniformity & smoothness
        uniformity(kk,ii+1)=sum(yBin.*yBin);
        variance=std(rawData{ii}(:))^2;
        smoothness(kk,ii+1)=1-(1/(1+variance)^2);
    end
end
close(h_wait)
uniformity(:,2:end)= ...
    IOS_time_gui_filter(uniformity(:,2:end),handles,'vertical');
handles.uniformity_meanOfgroup=uniformity;
smoothness(:,2:end)= ...
    IOS_time_gui_filter(smoothness(:,2:end),handles,'vertical');
handles.smoothness_meanOfgroup=smoothness;
handles.groupQuantity=groupQuantity;
guidata(hObject,handles)
% --------------------------------------------------------------------
function me_texture_show_Callback(hObject, eventdata, handles)
% hObject    handle to me_texture_show (see GCBO)
rawData=handles.smoothness_meanOfgroup;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(132);
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
for ii=1:length(pointPst)
    for jj=1:1
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
        legendName{vv-1}=['group',num2str(ii)];
    end
end
legend(legendName)
title('smoothness of group')
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;
handles.smoothness_meanOfgroup=rawData;
xlabel(timeCourseUnit)
hold off;
%%
rawData=handles.uniformity_meanOfgroup;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(133);
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
for ii=1:length(pointPst)
    for jj=1:1
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
        legendName{vv-1}=['group',num2str(ii)];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;
title('uniformity of group')
handles.uniformity_meanOfgroup=rawData;
xlabel(timeCourseUnit)
hold off;
guidata(hObject,handles)
% --------------------------------------------------------------------
function me_texture_export_Callback(hObject, eventdata, handles)
% hObject    handle to me_texture_export (see GCBO)
cd(handles.filePath)
result=get(handles.et_positiveNegativeFolder,'string');
resultFolder=fullfile(handles.filePath,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
[file,filePath]=uiputfile('uniformityAndSmoothness.mat');
cd(filePath)
file=file(1:end-4);
uniformity_meanOfgroup=handles.uniformity_meanOfgroup;
smoothness_meanOfgroup=handles.smoothness_meanOfgroup;
pointPst=handles.pointPst;
colorTmp=handles.colorTmp;
save([file,'_uniformity.txt'],'uniformity_meanOfgroup','-ASCII');
save([file,'_smoothness.txt'],'smoothness_meanOfgroup','-ASCII');
save([file,'.mat'],'uniformity_meanOfgroup','smoothness_meanOfgroup', ...
    'pointPst','colorTmp')
figure(132)
saveas(gcf,'smoothness','tiff')
figure(133)
saveas(gcf,'uniformity','tiff')
% --------------------------------------------------------------------
function me_texture_load_Callback(hObject, eventdata, handles)
% hObject    handle to me_texture_load (see GCBO)
[file,filePath]=uigetfile('uniformity*.mat');
load(fullfile(filePath,file))
handles.uniformity_meanOfgroup=uniformity_meanOfgroup;
handles.smoothness_meanOfgroup=smoothness_meanOfgroup;
handles.pointPst=pointPst;
handles.colorTmp=colorTmp;
guidata(hObject,handles)
function et_texture_edge_Callback(hObject, eventdata, handles)
% hObject    handle to et_texture_edge (see GCBO)
% et_texture_edge as text
%        str2double(get(hObject,'String')) returns contents of et_texture_edge as a double
function et_texture_edge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_texture_edge (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --------------------------------------------------------------------
function me_texture_singleHist_Callback(hObject, eventdata, handles)
% hObject    handle to me_texture_singleHist (see GCBO)
%%
Untitled_8_Callback(hObject, eventdata, handles)
% --- Executes on button press in cb_getIntensity.
function cb_getIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to cb_getIntensity (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of cb_getIntensity
% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% --------------------------------------------------------------------
function me_run_max_Callback(hObject, eventdata, handles)
% hObject    handle to me_run_max (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
rawData=zeros(p,sum(ROI_N(:))+1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
maxNumber=eval(get(handles.et_maxNumber,'string'));
h_wait=waitbar(0,'wait');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    img=differentTypeRead(file,fileType);
    vv=1;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
            vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            ROItmp=img(bw);
            dataTmp=sort(ROItmp(:),'descend');
            rawData(kk,vv)=mean(dataTmp(1:maxNumber));
        end
    end
end
close(h_wait)
rawData(:,2:end)= ...
    IOS_time_gui_filter(rawData(:,2:end),handles,'vertical');
handles.rawData_maxOfROI=rawData;
guidata(hObject,handles)
% --------------------------------------------------------------------
function me_resultShow_max_Callback(hObject, eventdata, handles)
% hObject    handle to me_resultShow_max (see GCBO)
rawData=handles.rawData_maxOfROI;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(103);
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj)];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;
handles.rawData_maxOfROI=rawData;
title('max')
xlabel(timeCourseUnit)
hold off;
guidata(hObject,handles)
% --------------------------------------------------------------------
function me_save_max_Callback(hObject, eventdata, handles)
% hObject    handle to me_save_max (see GCBO)
cd(handles.filePath)
result=get(handles.et_pearsonFolder,'string');
resultFolder=fullfile(handles.filePath,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
[file,filePath]=uiputfile('maxOfGroup.txt');
cd(filePath)
maxOfROI=handles.rawData_maxOfROI;
save(file,'maxOfROI','-ASCII');
% --------------------------------------------------------------------
function me_load_max_Callback(hObject, eventdata, handles)
% hObject    handle to me_load_max (see GCBO)
[file,filePath]=uigetfile('*.txt');
cd(filePath)
handles.rawData_maxOfROI=load(file);
guidata(hObject,handles)
% --------------------------------------------------------------------
function me_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to me_rotation (see GCBO)
% --------------------------------------------------------------------
function me_90_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to me_90_rotation (see GCBO)

file=handles.file;
filePath=handles.filePath;
cd(filePath)
fileType=handles.fileType;
result='rotate';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end
result=file;
% [result,resultPath]=uiputfile('*.*');
% cd(resultPath)
%differentTypeWrite(imgTmp{ii},file,fileType,handles.imgDepth);
% if get(handles.cb_rotation,'value')
    cd(handles.filePath)
    fileName=handles.fileName;
    p=length(fileName);
    h_wait=waitbar(0,'please wait');
%     if strcmp(handles.fileType,'mat')~=1
        for ii=1:p
            waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
            file=fileName(ii).name;
            img=differentTypeRead(fullfile(filePath,file),fileType);
            if size(img,3)==1
                differentTypeWrite(img.',fullfile(resultPath,file),fileType,handles.imgDepth);
            elseif size(img,3)==3
                img1=img(:,:,1).';img2=img(:,:,2).';img3=img(:,:,3).';
                img4=zeros(size(img1,1),size(img1,2),3);
                img4(:,:,1)=img1;
                img4(:,:,2)=img2;
                img4(:,:,3)=img3;
                differentTypeWrite(img4,fullfile(resultPath,file),fileType,handles.imgDepth);
            end
%                     img=imread(file);
        end
%     end
    close(h_wait)
% end
% --------------------------------------------------------------------
function me_upsideDown_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to me_upsideDown_rotation (see GCBO)
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
[result,resultPath]=uiputfile(file);
% cd(resultPath)
%differentTypeWrite(imgTmp{ii},file,fileType,handles.imgDepth);
% if get(handles.cb_rotation,'value')
    cd(handles.filePath)
    fileName=handles.fileName;
    p=length(fileName);
    h_wait=waitbar(0,'please wait');
%     if strcmp(handles.fileType,'mat')~=1
        for ii=1:p
            waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
            file=fileName(ii).name;
            img=differentTypeRead(fullfile(filePath,file),fileType);
            differentTypeWrite(img(end:-1:1,:,:),fullfile(resultPath,file),fileType,handles.imgDepth);
%                     img=imread(file);
        end
%     end
    close(h_wait)
% --------------------------------------------------------------------
function me_leftToRight_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to me_leftToRight_rotation (see GCBO)
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
[result,resultPath]=uiputfile(file);
% cd(resultPath)
%differentTypeWrite(imgTmp{ii},file,fileType,handles.imgDepth);
% if get(handles.cb_rotation,'value')
    cd(handles.filePath)
    fileName=handles.fileName;
    p=length(fileName);
    h_wait=waitbar(0,'please wait');
%     if strcmp(handles.fileType,'mat')~=1
        for ii=1:p
            waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
            file=fileName(ii).name;
            img=differentTypeRead(fullfile(filePath,file),fileType);
            differentTypeWrite(img(:,end:-1:1,:),fullfile(resultPath,file),fileType,handles.imgDepth);
%                     img=imread(file);
        end
%     end
    close(h_wait)
% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% --------------------------------------------------------------------
function me_zoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to me_zoomIn (see GCBO)
axes(handles.axes_rawImg) %#ok<*MAXES>
zoom on;
guidata(hObject,handles)
% --------------------------------------------------------------------
function me_zoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to me_zoomOut (see GCBO)
axes(handles.axes_rawImg)
zoom out;
zoom off;
guidata(hObject,handles)
% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
uiopen('*.txt')
function et_inline_texture_Callback(hObject, eventdata, handles)
% hObject    handle to et_inline_texture (see GCBO)
% et_inline_texture as text
%        str2double(get(hObject,'String')) returns contents of et_inline_texture as a double
function et_inline_texture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_inline_texture (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
handles.newIntensityGroup=inline(get(handles.et_inline_texture,'string'));
% get the bins.
edgeBin=eval(get(handles.et_texture_edge,'string'));
%pBin=length(edgeBin);
centerBin=edgeBin(1:end-1)+edgeBin(2:end);
centerBin=centerBin/2;
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
pointPst=handles.pointPst;
%p=length(fileName);
p=1;
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
yBin=zeros(length(centerBin),groupN);
ROI_N=zeros(groupN,1);
rawData=cell(groupN,1);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        bw=roipoly(img(:,:,1),cx,cy);
        ROI_N(ii)=ROI_N(ii)+sum(bw(:));
    end
%     minData=min();
    rawData{ii}=zeros(ROI_N(ii),1);
end
groupQuantity=ROI_N;
%h_wait=waitbar(0,'wait');
pp_trash=length(pointPst);
uniformity=zeros(p,pp_trash);
smoothness=zeros(p,pp_trash);
yBinSave=cell(pp_trash,1);
for kk=1:p
    %waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    %file=fileName(kk).name;
  %  img=differentTypeRead(file,fileType);
   % vv=1;
    for ii=1:length(pointPst)
        startIndex=1;
        for jj=1:length(pointPst{ii})
           % vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            ROItmp=img(bw);
            disp(['%%%%%%%%%%%%%%%%%', ...
                'min',num2str(min(ROItmp)),'max',num2str(max(ROItmp))])
            ROI_length=length(ROItmp);
            rawData{ii}(startIndex:startIndex+ROI_length-1)=ROItmp;
            startIndex=startIndex+ROI_length;
        end
        rawData{ii}=handles.newIntensityGroup(rawData{ii});
        % histogram
        yBinTmp=histc(rawData{ii}(:),edgeBin);
        %normalized
        yBinTmp=yBinTmp/groupQuantity(ii);
        % add the bin on the edge to the last
        yBinTmp(end-1)=yBinTmp(end-1)+yBinTmp(end);
        yBin=yBinTmp(1:end-1);
                uniformity(kk,ii)=sum(yBin.*yBin);
        variance=std(rawData{ii}(:))^2;
        smoothness(kk,ii)=1-(1/(1+variance)^2);
        yBinSave{ii}=yBin;
 %       percentTotal=sum(yBin(:))/mm/nn;
 %% uniformity & smoothness
    end
end
transferData.rawData=rawData;
transferData.centerBin=centerBin;
transferData.yBinSave=yBinSave;
transferData.uniformity=uniformity;
transferData.smoothness=smoothness;
handles.transferData=transferData;
handles.gui3_h=figure(131);close(handles.gui3_h);handles.gui3_h=figure(131);
set(handles.gui3_h,'unit','pixel','position',[31    94   648+100   464+100])
gui3String=cell(groupN,1);
% mean, max,min, smoothness, uniformity,group
handles.gui3_uitable1= uitable(handles.gui3_h, ...
    'data',cell(8,1), ...
    'position',[ 416+50   290   184+50   164+50], ...
    'rowName',{'group','uniformity','smoothness','max','mean','min','pixelQuantity','counted%'});
handles.gui3_axe1=axes('unit','pixel','position',[50   127   401   337]);
%set(handles.gui3_axe1,'parent',handles.gui3_h)
for ii=1:groupN
    gui3String{ii}=['group',num2str(ii)];
end
guidata(hObject,handles)
handles.gui3_listbox1 = uicontrol(handles.gui3_h,'Style','listbox', ...
    'string',gui3String, ...
    'position',[440+50 48 141 231], ...
    'callback',{@gui3_listbox1_callback,guidata(hObject)});
guidata(hObject,handles)
function gui3_listbox1_callback(hObject,eventdata,handles)
group=get(hObject,'Value');
transferData=handles.transferData;
rawData=transferData.rawData;
centerBin=transferData.centerBin;
yBin=transferData.yBinSave;
uniformity=transferData.uniformity;
smoothness=transferData.smoothness;
img=rawData{group};
img=double(img);
imgQuality=cell(8,1);
imgQuality{1}=num2str(group);
imgQuality{2}=num2str(uniformity(group));
imgQuality{3}=num2str(smoothness(group));
imgQuality{4}=num2str(max(img(:)));
imgQuality{5}=num2str(mean(img(:)));
imgQuality{6}=num2str(min(img(:)));
imgQuality{7}=num2str(length(img(:)));
imgQuality{8}=num2str(100*sum(yBin{group}));
set(handles.gui3_uitable1,'data',imgQuality)
plot(handles.gui3_axe1,centerBin,yBin{group})
% --- Executes on button press in cb_NoZero.
function cb_NoZero_Callback(hObject, eventdata, handles)
% hObject    handle to cb_NoZero (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of cb_NoZero
function et_step_Callback(hObject, eventdata, handles)
% hObject    handle to et_step (see GCBO)
% et_step as text
%        str2double(get(hObject,'String')) returns contents of et_step as a double
function et_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_step (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% --------------------------------------------------------------------
function me_stdOverTime_Callback(hObject, eventdata, handles)
% hObject    handle to me_stdOverTime (see GCBO)
filePath=handles.filePath;
file=handles.file;
fileName=handles.fileName;
p=length(fileName);
cd(filePath)
imgInfo=imfinfo(handles.fileName(1).name);
uint8Flag=imgInfo.BitDepth;
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1='STD';
resultPath1=fullfile(filePathTmp,resultName1);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
cd(filePath)
img=imread(file);
[m,n,nn_trash]=size(img);
imageStack=zeros(m,n,p,['uint',num2str(uint8Flag)]);
hh=waitbar(0,'please wait','name', ...
    'just wait');
for ii=1:p
  waitbar(ii/p,hh,[num2str(100*ii/p,'%03.1f'),'% completed, imread images']);
    img=imread(fileName(ii).name);
    imageStack(:,:,ii)=img(:,:,1);
end
close(hh)
hh=waitbar(0,'please wait','name', ...
    'just wait');
imgMean=zeros(m,n,'single');
imgStd=zeros(m,n,'single');
for ii=1:m
    waitbar(ii/p,hh,[num2str(100*ii/p,'%03.1f'),'% completed, calculating std']);
    imgMean(ii,:)=mean(single(imageStack(ii,:,:)),3);
    imgStd(ii,:)=std(single(imageStack(ii,:,:)),0,3);
end
close(hh)
cd(resultPath1)
[file1,filePath1]=uiputfile('stdImg.mat');
cd(filePath1)
save(file1,'imgStd');
imgMean(imgMean==0)=1;
imgStd_normalized=imgStd./imgMean;
file2=[file1(1:end-4),'_devidedByMean.mat'];
save(file2,'imgStd_normalized')
% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% --------------------------------------------------------------------
function me_saveCropImage_Callback(hObject, eventdata, handles)
% hObject    handle to me_saveCropImage (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
% imgInfo=imfinfo(handles.fileName(1).name);
% uint8Flag=imgInfo.BitDepth;
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1='cropped';
resultPath1=fullfile(filePathTmp,resultName1);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
cd(resultPath1)
[fileTrash,resultPath1]=uiputfile(handles.file);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
resultName=cell(groupN);
resultPath=cell(groupN);
cd(resultPath1)
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
    resultName{ii}=['group',num2str(ii)];
    resultPath{ii}=fullfile(resultPath1,resultName{ii});
    if exist(resultPath{ii})==7
    else
        mkdir(resultName{ii});
    end
end
xMin=zeros(groupN,1);
yMin=zeros(groupN,1);
xMax=zeros(groupN,1);
yMax=zeros(groupN,1);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        if jj==1
            xMinTmp=min(cx);
            yMinTmp=min(cy);
            xMaxTmp=max(cx);
            yMaxTmp=max(cy);
        end
        xMinTmp=min([xMinTmp,min(cx)]);
        yMinTmp=min([yMinTmp,min(cy)]);
        xMaxTmp=max([xMaxTmp,max(cx)]);
        yMaxTmp=max([yMaxTmp,max(cy)]);
    end
    xMin(ii)=round(xMinTmp);
    yMin(ii)=round(yMinTmp);
    xMax(ii)=round(xMaxTmp);
    yMax(ii)=round(yMaxTmp);
end
imgTmp=cell(groupN,1);
if strcmp(fileType,'mat')
    h_wait=waitbar(0,'wait');
    for kk=1:p
        waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
        file=fileName(kk).name;
        cd(filePath)
%         img=imread(file);
        img=differentTypeRead(file,fileType);
    %     vv=1;
        for ii=1:length(pointPst)
            imgTmp{ii}=img(yMin(ii):yMax(ii),xMin(ii):xMax(ii),:);
            cd(resultPath{ii})
            imgTmp2=imgTmp{ii};
            save(file,'imgTmp2')
     %       imwrite(imgTmp{ii},[file(1:end-length(fileType)),'.tif'],'tif')
     %       imwrite(imgTmp{ii},file,fileType)
     %       if strcmp('jpg',fileType)
      %          imwrite(imgTmp{ii},file,fileType);
     %       else
%             imwrite(imgTmp{ii},file,fileType,'compression','none');
     %       end
    %         for jj=1:length(pointPst{ii})
    %             vv=vv+1;
    %             cx=pointPst{ii}{jj}(:,1);
    %             cy=pointPst{ii}{jj}(:,2);
    %             bw=roipoly(img(:,:,1),cx,cy);
    %             ROItmp=img(bw);
    %             rawData(kk,vv)=mean(ROItmp(:));
    %         end
        end
    end
    close(h_wait)
else
    h_wait=waitbar(0,'wait');
    for kk=1:p
        waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
        file=fileName(kk).name;
        cd(filePath)
        img=imread(file);
    %     img=differentTypeRead(file,fileType);
    %     vv=1;
        for ii=1:length(pointPst)
            imgTmp{ii}=img(yMin(ii):yMax(ii),xMin(ii):xMax(ii),:);
            cd(resultPath{ii})
     %       imwrite(imgTmp{ii},[file(1:end-length(fileType)),'.tif'],'tif')
     %       imwrite(imgTmp{ii},file,fileType)
     %       if strcmp('jpg',fileType)
      %          imwrite(imgTmp{ii},file,fileType);
     %       else
            imwrite(imgTmp{ii},file,fileType);
     %       end
    %         for jj=1:length(pointPst{ii})
    %             vv=vv+1;
    %             cx=pointPst{ii}{jj}(:,1);
    %             cy=pointPst{ii}{jj}(:,2);
    %             bw=roipoly(img(:,:,1),cx,cy);
    %             ROItmp=img(bw);
    %             rawData(kk,vv)=mean(ROItmp(:));
    %         end
        end
    end
    close(h_wait)
end
pointPst=handles.pointPst;
colorTmp=handles.colorTmp;
cd ..
save('parameter.mat','pointPst','colorTmp')
% --------------------------------------------------------------------
function me_saveCropImageOnlyROI_Callback(hObject, eventdata, handles)
% hObject    handle to me_saveCropImageOnlyROI (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
% imgInfo=imfinfo(handles.fileName(1).name);
% uint8Flag=imgInfo.BitDepth;
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1='cropped_zeroPadding';
resultPath1=fullfile(filePathTmp,resultName1);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
resultName=cell(groupN);
resultPath=cell(groupN);
cd(resultPath1)
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
    resultName{ii}=['group',num2str(ii)];
    resultPath{ii}=fullfile(resultPath1,resultName{ii});
    if exist(resultPath{ii})==7
    else
        mkdir(resultName{ii});
    end
end
xMin=zeros(groupN,1);
yMin=zeros(groupN,1);
xMax=zeros(groupN,1);
yMax=zeros(groupN,1);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        if jj==1
            xMinTmp=min(cx);
            yMinTmp=min(cy);
            xMaxTmp=max(cx);
            yMaxTmp=max(cy);
        end
        xMinTmp=min([xMinTmp,min(cx)]);
        yMinTmp=min([yMinTmp,min(cy)]);
        xMaxTmp=max([xMaxTmp,max(cx)]);
        yMaxTmp=max([yMaxTmp,max(cy)]);
    end
    xMin(ii)=round(xMinTmp);
    yMin(ii)=round(yMinTmp);
    xMax(ii)=round(xMaxTmp);
    yMax(ii)=round(yMaxTmp);
end
imgTmp=cell(groupN,1);
h_wait=waitbar(0,'wait');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    cd(filePath)
    img=differentTypeRead(file,fileType);
    img2=img;
     img=differentTypeRead(file,fileType);
%     vv=1;
    for ii=1:length(pointPst)
        bwTmp=false(size(img));
        for jj=1:length(pointPst{ii})
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            bwTmp=or(bw,bwTmp);
        end
        img=img2;
        img(~bwTmp)=0;
        imgTmp{ii}=img(yMin(ii):yMax(ii),xMin(ii):xMax(ii),:);
        cd(resultPath{ii})
%         differentTypeWrite(img,file, fileType,imgDepth)
        differentTypeWrite(imgTmp{ii},file,fileType,handles.imgDepth);
%         imwrite(imgTmp{ii},file,fileType,'compression','none');
    end
end
close(h_wait)
function pushbutton85_Callback(hObject, eventdata, handles)
filePath=handles.filePath;
fileName=handles.fileName;
fileType=handles.fileType;
cd(filePath)
fileTmp=get(handles.et_fileName_file,'string');
[resultFile,resultPath]=uiputfile([fileTmp,'.',fileType]);
dotNO=strfind(resultFile,'.');
if isempty(dotNO)
    resultFile=[resultFile,'tif'];
    dotNO=strfind(resultFile,'.');
end
dotNO=dotNO(end);
file4=resultFile(dotNO+1:end);
file3='.';
file1=resultFile(1:dotNO-1);
if strcmp(resultPath,filePath)
       ButtonName = questdlg('are you sure want to replace raw images?', ...
           'choose', 'YES', 'NOT',  'NOT');
       if strcmp(ButtonName,'YES')
           h_wait=waitbar(0,'please wait');
          for ii=1:length(fileName)
              waitbar(ii/length(fileName),h_wait,[num2str(ii/length(fileName)*100,'%02d'),'%'])
              file2=num2str(ii,'%06d');
              file0=[file1,file2,file3,file4];
              file=fullfile(resultPath,file0);
              fileRaw=fullfile(filePath,fileName(ii).name);
              movefile(fileRaw,file);
          end
          close(h_wait)
       end
else
     h_wait=waitbar(0,'please wait');
    for ii=1:length(fileName)
        waitbar(ii/length(fileName),h_wait,[num2str(ii/length(fileName)*100,'%02d'),'%'])
        file2=num2str(ii,'%06d');
        file0=[file1,file2,file3,file4];
        file=fullfile(resultPath,file0);
        fileRaw=fullfile(filePath,fileName(ii).name);
        img=imread(fileRaw);
        imwrite(img,file,file4,'compression','none');
        %movefile(fileRaw,file);
    end
    close(h_wait)
end
function pushbutton86_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
%baselineNumber=eval(get(handles.et_dIOS_to,'string'));
binning=eval(get(handles.et_binning,'string'));
fileTypeLength=length(fileType);
selectedFileValue=get(handles.lb_curFileNames,'Value');
cd(filePath)
img2=differentTypeRead(file,fileType);
[mm,nn,nn_trash]=size(img2);
imgSingle=zeros(mm,nn,'single');
imgAvg=imgSingle;
cd(filePath)
 for jj=selectedFileValue:selectedFileValue+binning-1
    file=fileName(jj).name;
    img2=differentTypeRead(file,fileType);
if get(handles.cb_IOS_filter,'value')
    img2=differentTypeReadFilter_handles(img2,handles);
end    
    imgSingle=imgSingle+single(img2(:,:,1));
 end
 imgSingle=imgSingle/binning;
imgAvg=handles.baselineImg_IOS;
imgBack=handles.baselineImg_IOS_pre;
inverseNegative=get(handles.cb_inverseNegative,'value');
    if get(handles.cb_dF,'value')
        imgSave=(imgSingle-imgAvg);
    else
        imgSave=(imgSingle-imgAvg)./imgAvg;
    end
    
    %imgIOStmp=(imgAvg-baselineImg);
    if inverseNegative==1
        imgSave=abs(imgSave);
    end
    if get(handles.cb_biggerBack,'value')
        imgSave(eval(get(handles.et_IOSbigger,'string')))=0;
    end
    

% imgSave(handles.baselineImg_IOS_pre==0)=0;
convertFunction=get(handles.et_IOSinline,'string');
if isfield(handles,'ff_dIOS')
    if ishandle(handles.ff_dIOS)
        close(handles.ff_dIOS)
    end
end
handles.ff_dIOS=singleRetrieve_dIOS(imgSave,convertFunction);
guidata(hObject,handles)
function pushbutton87_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
%%
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1=get(handles.et_rawIOS_gray,'string');
resultName2=get(handles.et_rawIOS_pseudocolor,'string');

 if get(handles.cb_IOS_filter,'value')
    resultName1=[resultName1,handles.filterMethod2.saveName];
    resultName2=[resultName2,handles.filterMethod2.saveName];
 end
            
resultPath1=fullfile(filePathTmp,resultName1);
resultPath2=fullfile(filePathTmp,resultName2);
% cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultPath1);
end
if exist(resultPath2)==7
else
    mkdir(resultPath2);
end
cd(filePath)
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
fileName=handles.fileName;
p=length(fileName);
newIntensity=inline(get(handles.et_IOSinline,'string'));
contents = cellstr(get(handles.pp_IOScolormap,'String'));
colorSelectedTmp=contents(get(handles.pp_IOScolormap,'value'));
colorSelected=eval(colorSelectedTmp{1});
hh=waitbar(0,'please wait ...');
baselineImg=handles.baselineImg_IOS;
binningN=eval(get(handles.et_binning,'string'));
pp=floor(p/binningN);
img=differentTypeRead(fullfile(filePath,file),fileType);
if get(handles.cb_IOS_continuous,'value')
    imgAvg=zeros(size(img,1),size(img,2),binningN);
        for ii=1:binningN
            file=fileName(ii).name;
            img=differentTypeRead(fullfile(filePath,file),fileType);
             if get(handles.cb_IOS_filter,'value')
                img=differentTypeReadFilter_handles(img,handles);
            end           
            imgAvg(:,:,ii)=img;
        end 
        imgAvg(:,:,2:end)=imgAvg(:,:,1:end-1);
        for ii=binningN:p
                     waitbar(ii/p,hh, ...
                     [num2str(100*ii/p,'%03.1f'),'% completed']);           
            file=fileName(ii).name;
            img=differentTypeRead(fullfile(filePath,file),fileType);  
             if get(handles.cb_IOS_filter,'value')
                img=differentTypeReadFilter_handles(img,handles);
            end             
            imgAvg(:,:,1:end-1)=imgAvg(:,:,2:end);
            imgAvg(:,:,end)=img;
            imgIOStmp=(mean(imgAvg,3)-baselineImg)./baselineImg;
            img2=newIntensity(imgIOStmp);
            img3=uint8(img2);
            saveName1=[file(1:end-length(fileType)),'tif'];
            imwrite(img3,fullfile(resultPath1,saveName1),'tiff');
            saveName2=[file(1:end-length(fileType)),'png'];
            imgColor=ind2rgb(img3,colorSelected);
            imwrite(imgColor,fullfile(resultPath2,saveName2),'png')            
            
        end
else
    for ii=1:pp
                    waitbar(ii/pp,hh, ...
                     [num2str(100*ii/pp,'%03.1f'),'% completed']);
       imgAvg=img*0;
       cd(filePath)
        for jj=1:binningN
        file=fileName((ii-1)*binningN+jj).name;
        img=differentTypeRead(file,fileType);
             if get(handles.cb_IOS_filter,'value')
                img=differentTypeReadFilter_handles(img,handles);
            end         
        imgAvg=imgAvg+img;
        end
        imgAvg=imgAvg/binningN;
        imgIOStmp=(imgAvg-baselineImg)./baselineImg;
       % imgIOStmp=(imgAvg-baselineImg);
        img2=newIntensity(imgIOStmp);
        img3=uint8(img2);
%         cd(resultPath1)
        saveName1=[file(1:end-length(fileType)),'tif'];
        imwrite(img3,fullfile(resultPath1,saveName1),'tiff');
%         cd(resultPath2)
        saveName2=[file(1:end-length(fileType)),'png'];
        imgColor=ind2rgb(img3,colorSelected);
        imwrite(imgColor,fullfile(resultPath2,saveName2),'png')
        if get(handles.cb_movie,'value')
            M_gray(ii)=im2frame(img3,gray(256)); %#ok<*AGROW>
            M_color(ii)=im2frame(img3,colorSelected);
        end
    end
     if get(handles.cb_movie,'value')
        movie2avi(M_gray,fullfile(resultPath2,'movie'),'compression','none','fps',5);
        implay(M_gray)
        movie2avi(M_color,fullfile(resultPath2,'movie'),'compression','none','fps',5);
        % movie2avi(M2,fullfile(filePath,'movie2'),'compression','none','fps',5);
        implay(M_color)
     end    
end

close(hh)
guidata(hObject,handles)
% --- Executes on selection change in pp_horizontal.
function pp_horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to pp_horizontal (see GCBO)
% Hints: contents = cellstr(get(hObject,'String')) returns pp_horizontal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pp_horizontal
function pp_horizontal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_horizontal (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_mSequence_Callback(hObject, eventdata, handles)
% hObject    handle to et_mSequence (see GCBO)
% et_mSequence as text
%        str2double(get(hObject,'String')) returns contents of et_mSequence as a double
function et_mSequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_mSequence (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_mSequence.
function pb_mSequence_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mSequence (see GCBO)
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
selectLines=eval(get(handles.et_mSequence,'string'));
%%
hv=get(handles.pp_horizontal,'value');% 1: horizontal 2: vertical
cd(filePath)
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
fileName=handles.fileName;
p=length(fileName);
img=differentTypeRead(fullfile(filePath,file),fileType);
[mm,nn,n3]=size(img);
if hv==1
    sequenceImg=zeros(p,nn,'single');
%    file1='horizontal';
else
    sequenceImg=zeros(mm,p,'single');
 %   file1='vertical';
end
hh=waitbar(0,'please wait ...');
for ii=1:p
                waitbar(ii/p,hh, ...
                 [num2str(100*ii/p,'%03.1f'),'% completed,m_sequence']);
    file=fileName(ii).name;
    %cd(filePath)
    img=differentTypeRead(fullfile(filePath,file),fileType);
    img=differentTypeReadFilter_handles(img,handles);
    if hv==1 %horizontal
        sequenceImg(ii,:)=mean(img(selectLines,:,1),1);
    else % vertical
        sequenceImg(:,ii)=mean(img(:,selectLines,1),2);
    end
end
close(hh)
handles.sequenceImg=sequenceImg;
guidata(hObject,handles)
function pushbutton89_Callback(hObject, eventdata, handles)
%%
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1='m_sequence';
if get(handles.pp_spfilterMethod,'value')>1
    resultName1=[resultName1,handles.filterMethod2.saveName];
end
resultPath1=fullfile(filePathTmp,resultName1);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultPath1);
end
contents = cellstr(get(handles.pp_mSequence_color,'String'));
colorSelectedTmp=contents(get(handles.pp_mSequence_color,'value'));
colorSelected=eval(colorSelectedTmp{1});
binning=eval(get(handles.et_mSequence_binning,'string'));
hv=get(handles.pp_horizontal,'value');
%%
if hv==1
    file1='horizontal';
    [p,nn]=size(handles.sequenceImg);
    binP=floor(p/binning);
    sequenceImg2=zeros(binP,nn,'single');
    for ii=1:binP
        sequenceImg2(ii,:)=mean(handles.sequenceImg((ii-1)*binning+1:ii*binning,:),1);
    end
    sequenceImg2=IOS_time_gui_filter(sequenceImg2,handles,'vertical');
else
    file1='vertical';
    [mm,p]=size(handles.sequenceImg);
    binP=floor(p/binning);
    sequenceImg2=zeros(mm,binP,'single');
    for ii=1:binP
        sequenceImg2(:,ii)=mean(handles.sequenceImg(:,(ii-1)*binning+1:ii*binning),2);
    end
    sequenceImg3=sequenceImg2.';
    sequenceImg3=IOS_time_gui_filter(sequenceImg3,handles,'vertical');
    sequenceImg2=sequenceImg3.';
end
if hv==1
    if get(handles.cb_ms_dI,'value')
        pre=eval(get(handles.et_pre2,'string'));
        lineMean=mean(sequenceImg2(1:pre,:));
        lineMean(lineMean==0)=1;
        lineMean2=ones(p,1)*lineMean;
        sequenceImg2=(sequenceImg2-lineMean2)./lineMean2;
    else
    end
else
    if get(handles.cb_ms_dI,'value')
        pre=eval(get(handles.et_pre2,'string'));
        lineMean=mean(sequenceImg2(:,1:pre),2);
        lineMean(lineMean==0)=1;
        lineMean2=lineMean*ones(1,p);
        sequenceImg2=(sequenceImg2-lineMean2)./lineMean2;
    else
    end    
end
%newIntensity=inline(get(handles.et_mSequence_range,'string'));
% sequenceImg2=newIntensity(sequenceImg2);
%%
%newImgColor=ind2rgb(uint8(sequenceImg2),colorSelected);
mapRange=eval(get(handles.et_mSequence_range,'string'));
h110=figure(110);imshow(sequenceImg2,mapRange)
colormap(colorSelected)
colorbar
selectLines=eval(get(handles.et_mSequence,'string'));
file2=['_',num2str(selectLines(1)),'to',num2str(selectLines(end))];
sequenceImg=handles.sequenceImg;
cd(resultPath1)
[file0,filePath0]=uiputfile([file1,file2,'rawData.mat']);
cd(filePath0)
% saveas(h110,file0(1:end-4),'tiff')
save(file0,'sequenceImg');
img3=im2uint8(mat2gray(sequenceImg2,mapRange));
imgColor=ind2rgb(img3,colorSelected);
imwrite(imgColor,[file0(1:end-4),'noColorbar.tif'],'tif');
% if get(handles.pp_spfilterMethod,'value')>1
%     imwrite(im2uint8(mat2gray(sequenceImg2,mapRange)),colorSelected,[file0(1:end-4),'noColorbar',handles.filterMethod2.saveName,'.tif'],'tif');
% else
%     imwrite(im2uint8(mat2gray(sequenceImg2,mapRange)),colorSelected,[file0(1:end-4),'noColorbar.tif'],'tif');
% end

function et_mSequence_binning_Callback(hObject, eventdata, handles)
% hObject    handle to et_mSequence_binning (see GCBO)
% et_mSequence_binning as text
%        str2double(get(hObject,'String')) returns contents of et_mSequence_binning as a double
function et_mSequence_binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_mSequence_binning (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pp_mSequence_color.
function pp_mSequence_color_Callback(hObject, eventdata, handles)
% hObject    handle to pp_mSequence_color (see GCBO)
% Hints: contents = cellstr(get(hObject,'String')) returns pp_mSequence_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pp_mSequence_color
function pp_mSequence_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_mSequence_color (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_mSequence_range_Callback(hObject, eventdata, handles)
% hObject    handle to et_mSequence_range (see GCBO)
% et_mSequence_range as text
%        str2double(get(hObject,'String')) returns contents of et_mSequence_range as a double
function et_mSequence_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_mSequence_range (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton90_Callback(hObject, eventdata, handles)
contents = cellstr(get(handles.pp_mSequence_color,'String'));
colorSelectedTmp=contents(get(handles.pp_mSequence_color,'value'));
colorSelected=eval(colorSelectedTmp{1});
binning=eval(get(handles.et_mSequence_binning,'string'));
hv=get(handles.pp_horizontal,'value');
%%
if hv==1
    [p,nn]=size(handles.sequenceImg);
    binP=floor(p/binning);
    sequenceImg2=zeros(binP,nn,'single');
    for ii=1:binP
        sequenceImg2(ii,:)=mean(handles.sequenceImg((ii-1)*binning+1:ii*binning,:),1);
    end
    sequenceImg2=IOS_time_gui_filter(sequenceImg2,handles,'vertical');
else
    [mm,p]=size(handles.sequenceImg);
    binP=floor(p/binning);
    sequenceImg2=zeros(mm,binP,'single');
    for ii=1:binP
        sequenceImg2(:,ii)=mean(handles.sequenceImg(:,(ii-1)*binning+1:ii*binning),2);
    end
    sequenceImg3=sequenceImg2.';
    sequenceImg3=IOS_time_gui_filter(sequenceImg3,handles,'vertical');
    sequenceImg2=sequenceImg3.';
end
if hv==1
    if get(handles.cb_ms_dI,'value')
        pre=eval(get(handles.et_pre2,'string'));
        lineMean=mean(sequenceImg2(1:pre,:));
        lineMean(lineMean==0)=1;
        lineMean2=ones(p,1)*lineMean;
        sequenceImg2=(sequenceImg2-lineMean2)./lineMean2;
    else
    end
else
    if get(handles.cb_ms_dI,'value')
        pre=eval(get(handles.et_pre2,'string'));
        lineMean=mean(sequenceImg2(:,1:pre),2);
        lineMean(lineMean==0)=1;
        lineMean2=lineMean*ones(1,p);
        sequenceImg2=(sequenceImg2-lineMean2)./lineMean2;
    else
    end    
end
%newIntensity=inline(get(handles.et_mSequence_range,'string'));
% sequenceImg2=newIntensity(sequenceImg2);
%%
%newImgColor=ind2rgb(uint8(sequenceImg2),colorSelected);
mapRange=eval(get(handles.et_mSequence_range,'string'));
figure(110);imshow(sequenceImg2,mapRange)
colormap(colorSelected)
colorbar
function pushbutton91_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.mat');
load(fullfile(filePath,file))
handles.sequenceImg=sequenceImg;
guidata(hObject,handles)
function pushbutton92_Callback(hObject, eventdata, handles)
handles=IOS_time_parameterRead(handles);
pre1=handles.pre1;
pre2=handles.pre2;
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
fileName=handles.fileName;
cd(filePath);
img=differentTypeRead(file,fileType);
p=handles.p;
inverseNegative=get(handles.cb_inverseNegative,'value');
%% prestimulus mean
baselineImg=zeros(size(img));
h_wait=waitbar(0,'please wait');
for ii=pre1:pre2
    waitbar(ii/pre2,h_wait,[num2str(100*ii/pre2,'%04.1f'),'% complated: import baseline images'])
    file=fileName(ii).name;
    img=differentTypeRead(file,fileType);
    baselineImg=baselineImg+img;
end
baselineImg=baselineImg/(pre2-pre1+1);
if get(handles.cb_IOS_filter,'value')
    baselineImg=differentTypeReadFilter_handles(baselineImg,handles);
end
handles.baselineImg_IOS_pre=baselineImg;
baselineImg(baselineImg==0)=1;
% baselineImg(baselineImg<200)=0;
handles.baselineImg_IOS=baselineImg;
close(h_wait)
guidata(hObject,handles)
function et_localSD_x_Callback(hObject, eventdata, handles)
% hObject    handle to et_localSD_x (see GCBO)
% et_localSD_x as text
%        str2double(get(hObject,'String')) returns contents of et_localSD_x as a double
function et_localSD_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_localSD_x (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_localSD_y_Callback(hObject, eventdata, handles)
% hObject    handle to et_localSD_y (see GCBO)
% et_localSD_y as text
%        str2double(get(hObject,'String')) returns contents of et_localSD_y as a double
function et_localSD_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_localSD_y (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_localSD_name_Callback(hObject, eventdata, handles)
% hObject    handle to et_localSD_name (see GCBO)
% et_localSD_name as text
%        str2double(get(hObject,'String')) returns contents of et_localSD_name as a double
function et_localSD_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_localSD_name (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton93_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
fileTypeLength=length(fileType);
windowX=eval(get(handles.et_localSD_x,'string'));
windowY=eval(get(handles.et_localSD_y,'string'));
conversion=get(handles.et_localSD_measure,'string');
sizeY=round((windowY-1)/2)*2+1;
sizeX=round((windowX-1)/2)*2+1;
sizeYhalf=round((windowY-1)/2);
sizeXhalf=round((windowX-1)/2);
IOSresultPath=get(handles.et_IOSresultPath,'string');
result=get(handles.et_localSD_name,'string');
resultPath=fullfile(IOSresultPath,result);
cd(IOSresultPath)
if exist(resultPath)==7
else
    mkdir(result)
end
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[mm,nn,nn_trash]=size(img);
y1=sizeYhalf+1;y2=mm-sizeYhalf;
x1=sizeXhalf+1;x2=nn-sizeXhalf;
% roiY=[y1,1,mm,mm,y2,y2];
% roiX=[x1,1,1,nn,x2,x1];
% roiY2=[y1,1,1,mm,y2,y1];
% roiX2=[x1,1,nn,nn,x2,x2];
h_wait=waitbar(0,'completed');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(kk/p*100,'%03.1f'),'%completed'])
    file=fileName(kk).name;
    cd(filePath)
    img=differentTypeRead(file,fileType);
    %% core
    stdImg=zeros(mm,nn);
    for ii=round((windowY-1)/2)+1:mm-round((windowY-1)/2)
        for jj=round((windowX-1)/2)+1:nn-round((windowX-1)/2)
            windowImg=img(ii-sizeYhalf:ii+sizeYhalf,jj-sizeXhalf:jj+sizeXhalf);
            %% get localized SD
            a=windowImg(:);
            eval(conversion)
            stdImg(ii,jj)=b;
            %b=std(a)/mean(a);
           % stdImg(ii,jj)=std(windowImg(:))/mean(windowImg(:));
        end
        %% duplicating padding over margin area
%         stdImg=roifill(stdImg,roiX,roiY);
%         stdImg=roifill(stdImg,roiX2,roiY2);
        stdImg(:,1:x1)=stdImg(:,x1)*ones(1,x1);
        stdImg(:,x2:end)=stdImg(:,x2)*ones(1,x1);
        stdImg(1:y1,:)=ones(y1,1)*stdImg(y1,:);
        stdImg(y2:end,:)=ones(y1,1)*stdImg(y2,:);
    end
    cd(resultPath)
        saveName=[file(1:end-fileTypeLength),'mat'];
    save(saveName,'stdImg');
end
close(h_wait)
function et_localSD_measure_Callback(hObject, eventdata, handles)
% hObject    handle to et_localSD_measure (see GCBO)
% et_localSD_measure as text
%        str2double(get(hObject,'String')) returns contents of et_localSD_measure as a double
function et_localSD_measure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_localSD_measure (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_binningTime_Callback(hObject, eventdata, handles)
% hObject    handle to et_binningTime (see GCBO)
% et_binningTime as text
%        str2double(get(hObject,'String')) returns contents of et_binningTime as a double
function et_binningTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_binningTime (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function handles=pb_binning_Callback(hObject, eventdata, handles)
fileName=handles.fileName;
filePath=handles.filePath;
file=handles.file;
p=length(fileName);
fileType=handles.fileType;
fileTypeLength=length(fileType);
n=eval(get(handles.et_binningTime,'string'));
IOSresultPath=get(handles.et_IOSresultPath,'string');
%% file name for result
contents = cellstr(get(handles.pp_binningMeanMax,'String'));
name=['Bin_',num2str(n),contents{get(handles.pp_binningMeanMax,'value')}];
result0='BinOverT';
IOSresultPath=get(handles.et_IOSresultPath,'string');
resultPath0=fullfile(IOSresultPath,result0);
cd(IOSresultPath)
if exist(resultPath0)==7
else
    mkdir(result0)
end
cd(resultPath0)
%% subdirectory
result=name;
resultPath=fullfile(resultPath0,result);
cd(resultPath0)
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
[fileTmp,resultPath]=uiputfile(file,'note changing name does not work');
%%
t=floor(p/n);
nameGroup=cell(t,1);
for ii=1:t
    nameGroup{ii}=fileName(n*(ii-1)+1:n*(ii-1)+n);
end
if mod(p,n)~=0
    ii=t;
    nameGroup{ii}=fileName(n*(ii-1)+1:end);
end
cd(filePath)
img=differentTypeRead(file,fileType);
[mm,nn,nn3]=size(img);
hh=waitbar(0,'please wait','name', ...
    'binninig');
binningRange=eval(get(handles.et_binningRange,'string'));
binningMethod=get(handles.pp_binningMeanMax,'value');
% maxQuantity=eval(get(handles.et_maxNumber,'string'));
cd(filePath)
if strcmp(fileType,'mat')==0
imgInfo=imfinfo(fileName(1).name);
uint8Flag=imgInfo.BitDepth;
end
if binningMethod<=3
     for ii=1:t
             waitbar(ii/t,hh, ...
             [num2str(100*ii/t,'%03.1f'),'% completed'], ...
             'binning');
        cd(filePath)
        fileNameTmp=nameGroup{ii};
        pp=length(fileNameTmp);
        imgStack=zeros(mm,nn,pp,'single');
         for jj=1:pp
            img=differentTypeRead(fileNameTmp(jj).name,fileType);
            imgStack(:,:,jj)=img(:,:,1);
         end
         if length(binningRange)<n
             imgStack=sort(imgStack,3,'descend');
         elseif binningMethod==2 && maxQuantity>1
             imgStack=sort(imgStack,3,'descend');
         end
        if binningMethod==1 %% mean
            imgAvg=mean(imgStack(:,:,binningRange),3);
        elseif binningMethod==2 %% max
            if maxQuantity==1 && length(binningRange)==n
                imgAvg=max(imgStack(:,:,:),[],3);
            else
                imgAvg=mean(imgStack(:,:,binningRange(1):binningRange(1)+maxQuantity-1),3);
            end
        elseif binningMethod==3 %% median
            imgAvg=median(imgStack(:,:,binningRange),3);
        end
        %% imgStack=zeros(size(img,1),size(img,2),pp,'single');
        cd(resultPath)
        differentTypeWrite(imgAvg,fileNameTmp(end).name,fileType,handles.imgDepth);
        %%
    %     if strcmp(fileType,'mat')==0
    %
    %         imwrite(eval(['uint',num2str(uint8Flag),'(imgAvg)']),fileNameTmp(end).name,fileType,'compression','none')
    %     %%
    %     else
    %         saveName=[fileNameTmp(end).name(1:end-fileTypeLength),'mat'];
    %         save(saveName,'imgAvg');
    %     end
    end
    close(hh)
elseif binningMethod==4
    matlabpool
    if mod(p,n)~=0
        t=t+1;
    end
    parfor ii=1:n
%               waitbar(ii/n,hh, ...
%              [num2str(100*ii/n,'%03.1f'),'% completed'], ...
%              'binning');
%         cd(filePath)
%         img=differentTypeRead(fileName(ii).name,fileType);
        imgAvg=zeros(mm,nn,nn3);
        jj_All=0;
        cd(filePath)
        for jjTmp=1:t
            jj=ii+(jjTmp-1)*n;
            if jj<=p
            img=differentTypeRead(fileName(jj).name,fileType);
            imgAvg=imgAvg+double(img);
            jj_All=jj_All+1;
            end
%             differentTypeWrite(imgAvg,fileName(ii).name,fileType,handles.imgDepth);
        end
        if jj_All>0
            imgAvg=imgAvg/jj_All;
        end
        cd(resultPath)
        differentTypeWrite(imgAvg,fileName(ii).name(1:end-length(fileType)-1),'mat',handles.imgDepth)
    end
    matlabpool close
    close(hh)
end
%% max pp_binningMeanMax
function et_xBinning_Callback(hObject, eventdata, handles)
% hObject    handle to et_xBinning (see GCBO)
% et_xBinning as text
%        str2double(get(hObject,'String')) returns contents of et_xBinning as a double
function et_xBinning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_xBinning (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_yBinning_Callback(hObject, eventdata, handles)
% hObject    handle to et_yBinning (see GCBO)
% et_yBinning as text
%        str2double(get(hObject,'String')) returns contents of et_yBinning as a double
function et_yBinning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_yBinning (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_xBinningStart_Callback(hObject, eventdata, handles)
% hObject    handle to et_xBinningStart (see GCBO)
% et_xBinningStart as text
%        str2double(get(hObject,'String')) returns contents of et_xBinningStart as a double
function et_xBinningStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_xBinningStart (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_yBinningStart_Callback(hObject, eventdata, handles)
% hObject    handle to et_yBinningStart (see GCBO)
% et_yBinningStart as text
%        str2double(get(hObject,'String')) returns contents of et_yBinningStart as a double
function et_yBinningStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_yBinningStart (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton95_Callback(hObject, eventdata, handles)
fileName=handles.fileName;
filePath=handles.filePath;
file=handles.file;
p=length(fileName);
fileType=handles.fileType;
fileTypeLength=length(fileType);
%% binning size
xBinning=eval(get(handles.et_xBinning,'string'));
yBinning=eval(get(handles.et_yBinning,'string'));
%% frame rate unit HZ
frameRate=1000;
%% ROI
xStart=eval(get(handles.et_xBinningStart,'string'));
xEnd=eval(get(handles.et_xBinningEnd,'string'));
yStart=eval(get(handles.et_yBinningStart,'string'));
yEnd=eval(get(handles.et_yBinningEnd,'string'));
%% verify whether ROI fits
if mod(yEnd-yStart+1,yBinning)~=0
    yEnd=floor((yEnd-yStart+1)/yBinning)*yBinning-1+yStart;
end
if mod(xEnd-xStart+1,xBinning)~=0
    xEnd=floor((xEnd-xStart+1)/xBinning)*xBinning-1+xStart;
end
new_xsize=(xEnd-xStart+1)/xBinning; % image size after pixel binning
new_ysize=(yEnd-yStart+1)/yBinning; % image size after pixel binning
%% file name for result
name=['x',num2str(xBinning),'y',num2str(yBinning), ...
    'ROIx',num2str(xStart),'to',num2str(xEnd),'ROIy',num2str(yStart),'to',...
    num2str(yEnd)];
result0='BinSpatl';
IOSresultPath=get(handles.et_IOSresultPath,'string');
resultPath0=fullfile(IOSresultPath,result0);
cd(IOSresultPath)
if exist(resultPath0)==7
else
    mkdir(result0)
end
cd(resultPath0)
%% subdirectory
result=name;
resultPath=fullfile(resultPath0,result);
cd(resultPath0)
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
[fileTmp,resultPath]=uiputfile(file,'note changing name does not work');
%%
h_wait=waitbar(0,'please wait');
p_new=p;
cd(filePath)
if strcmp(fileType,'mat')==0
    imgInfo=imfinfo(fileName(1).name);
    uint8Flag=imgInfo.BitDepth;
end
for ii=1:p_new
    waitbar(ii/p_new,h_wait,[num2str(100*ii/p_new,'%04.1f'),'% completed'])
    cd(filePath)
    img2=differentTypeRead(fileName(ii).name,fileType);
    img=image_binning(img2,xBinning,yBinning,xStart,yStart,new_xsize,new_ysize);
    cd(resultPath)
    if strcmp(fileType,'mat')==0
        imwrite(eval(['uint',num2str(uint8Flag),'(img)']),fileName(ii).name,fileType,'compression','none')
    else
      saveName=[fileName(ii).name(1:end-fileTypeLength),'mat'];
      save(saveName,'img');
    end
end
close(h_wait)
function imgBinning=image_binning(img,xBinning,yBinning,xStart,yStart,new_xsize,new_ysize)
imgBinning=zeros(new_ysize,new_xsize,'single');
if xBinning==1 && yBinning<=10
    %% speed up the processing for line-modulaiton OCT
    img0=img(yStart+(1-1)*yBinning:yStart+new_ysize*yBinning-1, ...
                xStart+(1-1)*xBinning:xStart+new_xsize*xBinning-1);
    img1=img0;
    img2=img0;
    for ii=2:yBinning
        img2(ii:end,:)=img0(1:end-ii+1,:);
        img1=img1+img2;
    end
    imgBinning=img1(1:yBinning:end,:)/yBinning;
elseif yBinning==1 && xBinning<=10
    img0=img(yStart+(1-1)*yBinning:yStart+new_ysize*yBinning-1, ...
                xStart+(1-1)*xBinning:xStart+new_xsize*xBinning-1);
    img1=img0;
    img2=img0;
    for ii=2:xBinning
        img2(:,ii:end)=img0(:,1:end-ii+1);
        img1=img1+img2;
    end
    imgBinning=img1(:,1:xBinning:end)/xBinning;
else
    for ii=1:new_ysize
        for jj=1:new_xsize
            imgROI=img(yStart+(ii-1)*yBinning:yStart+ii*yBinning-1, ...
                xStart+(jj-1)*xBinning:xStart+jj*xBinning-1);
            imgBinning(ii,jj)=mean(imgROI(:));
        end
    end
end
%function et_xBinningEnd_Callback(hObject, eventdata, handles)
% hObject    handle to et_xBinningEnd (see GCBO)
% et_xBinningEnd as text
%        str2double(get(hObject,'String')) returns contents of et_xBinningEnd as a double
% function et_xBinningEnd_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to et_xBinningEnd (see GCBO)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function et_yBinningEnd_Callback(hObject, eventdata, handles)
% % hObject    handle to et_yBinningEnd (see GCBO)
% % et_yBinningEnd as text
% %        str2double(get(hObject,'String')) returns contents of et_yBinningEnd as a double
% function et_yBinningEnd_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to et_yBinningEnd (see GCBO)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
function pushbutton96_Callback(hObject, eventdata, handles)
% --- Executes on selection change in popupmenu15.
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pp_binningMeanMax.
function pp_binningMeanMax_Callback(hObject, eventdata, handles)
% hObject    handle to pp_binningMeanMax (see GCBO)
% Hints: contents = cellstr(get(hObject,'String')) returns pp_binningMeanMax contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pp_binningMeanMax
function pp_binningMeanMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_binningMeanMax (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit71_Callback(hObject, eventdata, handles)
% hObject    handle to edit71 (see GCBO)
% edit71 as text
%        str2double(get(hObject,'String')) returns contents of edit71 as a double
function edit71_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit71 (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit72_Callback(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% edit72 as text
%        str2double(get(hObject,'String')) returns contents of edit72 as a double
function edit72_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pp_spfilterMethod.
function handles=pp_spfilterMethod_Callback(hObject, eventdata, handles)
% hObject    handle to pp_spfilterMethod (see GCBO)
handles=pp_2filterMethod_outside(hObject,handles);
if handles.filterMethod2.name==5
    filterMethod2=handles.filterMethod2;
    if isfield(filterMethod2,'answer')
        def=filterMethod2.answer;
    else
        def = {'4','0.5'};
        
    end
    
%     set(handles.text94,'string','order')
%     set(handles.p_2nonFilter,'visible','off');
    prompt = {'fitting order:','fcut'};
    dlg_title = 'Input for butterworth filter on Foureir Domain';
    num_lines = 1;
    
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    handles.filterMethod2.order=eval(answer{1});
    handles.filterMethod2.fcut=eval(answer{2});
    handles.filterMethod2.answer=answer;
%     handles.filterMethod.detrend_base=answer{3};s
    handles.filterMethod2.saveName=['butter_order',answer{1},'_fcut',answer{2}];
end

guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns pp_spfilterMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pp_spfilterMethod
function handles=pp_2filterMethod_outside(hObject,handles)
handles.filterMethod2.name=get(handles.pp_spfilterMethod,'Value')-1;
handles=filterInfoRead2_spatial(handles);
nameTmp=handles.filterMethod2.name;
if nameTmp==0 %non filter
    set(handles.p_2nonFilter,'visible','off');
    handles.filterMethod2.saveName=['noneFilter'];
elseif nameTmp==1% gaussian filter
    set(handles.p_2nonFilter,'visible','on');
    handles.filterMethod2.saveName=['Gauss_siz',get(handles.et_2sizNo,'string'), ...
        'sigma',get(handles.et_2sigmaNo,'string')];
    set(handles.p_2wiener,'visible','on')
elseif nameTmp==2 % median filter
   set(handles.p_2nonFilter,'visible','on');
   set(handles.p_2wiener,'visible','off');
    handles.filterMethod2.saveName=['Median_siz',get(handles.et_2sizNo,'string')];
elseif nameTmp==3
    set(handles.p_2nonFilter,'visible','on');
    set(handles.p_2wiener,'visible','off');
    handles.filterMethod2.saveName=['Wiener_siz',get(handles.et_2sizNo,'string')];
elseif nameTmp==4
    set(handles.p_2nonFilter,'visible','on');
    set(handles.p_2wiener,'visible','off');
    handles.filterMethod2.saveName=['max_siz',get(handles.et_2sizNo,'string')];
elseif nameTmp==6

    
end
function handles=filterInfoRead2_spatial(handles)
handles.filterMethod2.sizNo=str2num(get(handles.et_2sizNo,'string')); %#ok<*ST2NM>
%handles.sizSuper=str2num(get(handles.et_sizSuper,'string'));
handles.filterMethod2.sigmaNo=str2num(get(handles.et_2sigmaNo,'string'));
%handles.sigmaSuper=str2num(get(handles.et_sigmaSuper,'string'));
handles.filterMethod2.name=get(handles.pp_spfilterMethod,'Value')-1;
function pp_spfilterMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_spfilterMethod (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_2sizNo_Callback(hObject, eventdata, handles)
% hObject    handle to et_2sizNo (see GCBO)
% et_2sizNo as text
%        str2double(get(hObject,'String')) returns contents of et_2sizNo as a double
function et_2sizNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_2sizNo (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_2sigmaNo_Callback(hObject, eventdata, handles)
% hObject    handle to et_2sigmaNo (see GCBO)
% et_2sigmaNo as text
%        str2double(get(hObject,'String')) returns contents of et_2sigmaNo as a double
function et_2sigmaNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_2sigmaNo (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function img=differentTypeReadFilter_handles(img0,handles)
if get(handles.pp_spfilterMethod,'value')==1
    img=img0;
else
%     handles=filterInfoRead2_spatial(handles);
    filterMethod2=handles.filterMethod2;
    img=differentTypeReadFilter(img0,filterMethod2);    
end

function handles=pushbutton97_Callback(hObject, eventdata, handles)
fileName=handles.fileName;
filePath=handles.filePath;
file=handles.file;
p=length(fileName);
fileType=handles.fileType;
fileTypeLength=length(fileType);
%% file name for result
handles=pp_2filterMethod_outside(hObject,handles);
handles=filterInfoRead2_spatial(handles);
name=handles.filterMethod2.saveName;
result0='filter';
IOSresultPath=get(handles.et_IOSresultPath,'string');
resultPath0=fullfile(IOSresultPath,result0);
cd(IOSresultPath)
if exist(resultPath0)==7
else
    mkdir(result0)
end
cd(resultPath0)
%% subdirectory
result=name;
resultPath=fullfile(resultPath0,result);
cd(resultPath0)
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
[fileTmp,resultPath]=uiputfile(file,'note changing name does not work');
handles=filterInfoRead2_spatial(handles);%img=imgFilter(img,filterMethod2)
filterMethod2=handles.filterMethod2;
sizNo=filterMethod2.sizNo;
sigmaNo=filterMethod2.sigmaNo;
nameNo=filterMethod2.name;
%%
h_wait=waitbar(0,'please wait');
p_new=p;
cd(filePath)
if strcmp(fileType,'mat')==0
    imgInfo=imfinfo(fileName(1).name);
    uint8Flag=imgInfo.BitDepth;
end
for ii=1:p_new
    waitbar(ii/p_new,h_wait,[num2str(100*ii/p_new,'%04.1f'),'% completed'])
    img0=differentTypeRead(fullfile(filePath,fileName(ii).name),fileType);
    img=differentTypeReadFilter(img0,filterMethod2);
    
    differentTypeWrite(img,fullfile(resultPath,fileName(ii).name), fileType,handles.imgDepth)
%     if strcmp(fileType,'mat')==0
%         imwrite(eval(['uint',num2str(uint8Flag),'(img)']),fileName(ii).name,fileType,'compression','none')
%     else
%       saveName=[fileName(ii).name(1:end-fileTypeLength),'mat'];
%       save(saveName,'img');
%     end
end
close(h_wait)
cd(resultPath)
function et_binningRange_Callback(hObject, eventdata, handles)
% hObject    handle to et_binningRange (see GCBO)
% et_binningRange as text
%        str2double(get(hObject,'String')) returns contents of et_binningRange as a double
function et_binningRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_binningRange (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rePath_moveCreate
pathWhich=which('IOS_Software');
dotNO=strfind(pathWhich,filesep);
dotNO=dotNO(end);
pathWhich=pathWhich(1:dotNO-1);
filePath=pathWhich;
file='movCreate';
filePath=fullfile(filePath,file);
path(filePath,path);
function pb_saveMovie_Callback(hObject, eventdata, handles)
% rePath_moveCreate;
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
subPath = split(filePath,filesep);
subPath = subPath(~cellfun('isempty',subPath));
if ~exist(handles.filePath_movie,'dir')
    mkdir(handles.filePath_movie);
end

%%
% filePathTmp=get(handles.et_IOSresultPath,'string');
% resultName1='movie';
% resultPath1=fullfile(filePathTmp,resultName1);
% cd(filePathTmp)
% if exist(resultPath1)==7
% else
%     mkdir(resultName1);
% end
% cd(filePath)
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
fileName=handles.fileName;
p=length(fileName);
movieFPS=eval(get(handles.et_fps,'string'));
newIntensity=inline(get(handles.et_inline,'string'));
contents = cellstr(get(handles.pm_colormap_raw,'String'));
colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
colorSelected=eval(colorSelectedTmp{1});

 % cd(filePath)
%  hh=waitbar(0,'please wait ...');
% for ii=1:p
%                 waitbar(ii/p,hh, ...
%                  [num2str(100*ii/p,'%03.1f'),'% completed']);
%     file=fileName(ii).name;
%     img=differentTypeRead(fullfile(filePath,file),fileType);
%     img2=newIntensity(img);
%     if get(handles.cb_eval,'value')
%         x=img2;
%         eval(get(handles.et_eval,'string')) ;
%         img2=x;
%     end
%     img3=uint8(img2);
%         M_color(ii)=im2frame(img3,colorSelected);
% end
% close(hh)
%     cd(resultPath1)
%     movieName='movie.avi';
% % [movieName,resultPath1]=uiputfile('movie.avi');
% movieName=movieName(1:end-4);
% cd(resultPath1)
% 
%     movie2avi(M_color,fullfile(resultPath1,[movieName,'2.avi']),'compression','none','fps',movieFPS);

    % movie2avi(M2,fullfile(filePath,'movie2'),'compression','none','fps',5);
%     implay(M_color,movieFPS)
%%
movieName=[subPath{end} '-movie.mp4'];
[movieName,resultPath1]=uiputfile(fullfile(handles.filePath_movie,movieName));
vidObj = VideoWriter([resultPath1,filesep,movieName],'MPEG-4');
%     vidObj = VideoWriter([movieName(1:end-4),'_2'],'mpeg-4');
vidObj.FrameRate=movieFPS;
vidObj.Quality=100;
%     vidObj.LosslessCompression=1;
open(vidObj);
hh=waitbar(0,'please wait ...','Name','Saving the movie');
% cd(filePath)
for ii=1:p
                waitbar(ii/p,hh, ...
                 [num2str(100*ii/p,'%03.1f'),'% completed']);
    file=fileName(ii).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    img2=newIntensity(img);
    if get(handles.cb_eval,'value')
        x=img2;
        eval(get(handles.et_eval,'string')) ;
        img2=x;
    end
    if size(img2,3)==3

        currFrame=im2frame(uint8(img2));
               writeVideo(vidObj,currFrame);             
    else
    img3=uint8(img2(:,:));
        currFrame=im2frame(img3,colorSelected);
               writeVideo(vidObj,currFrame);        
    end

end
close(hh)
close(vidObj);
% %% mov
% movieName='moviemov.mov';
% % cd(resultPath1)
% % [movieName,resultPath1]=uiputfile(movieName);
% cd(resultPath1)
% vidObj  = QTWriter(movieName,'MovieFormat','Photo JPEG','Quality',100);
% % vidObj  = QTWriter(movieName,'MovieFormat','Photo TIFF');
% % vidObj  = QTWriter(movieName);
% %     vidObj = VideoWriter([movieName(1:end-4),'_2'],'mpeg-4');
% vidObj.FrameRate=movieFPS;
% % vidObj.Quality=100;
% %     vidObj.LosslessCompression=1;
% % open(vidObj);
% hh=waitbar(0,'please wait ...');
% cd(filePath)
% for ii=1:p
%                 waitbar(ii/p,hh, ...
%                  [num2str(100*ii/p,'%03.1f'),'% completed']);
%     file=fileName(ii).name;
%     img=differentTypeRead(file,fileType);
%     img2=newIntensity(img);
%     if get(handles.cb_eval,'value')
%         x=img2;
%         eval(get(handles.et_eval,'string')) ;
%         img2=x;
%     end
%     if size(img2,3)==3
% 
%         currFrame=im2frame(uint8(img2));
%                writeMovie(vidObj,currFrame);             
%     else
%     img3=uint8(img2(:,:));
%         currFrame=im2frame(img3,colorSelected);
%                writeMovie(vidObj,currFrame);        
%     end
% 
% end
% close(hh)
% close(vidObj);
% 
% %%
% movieName='movieAvi.avi';
% cd(resultPath1)
% % [movieName,resultPath1]=uiputfile(movieName);
% % cd(resultPath1)
% vidObj = VideoWriter(movieName);
% %     vidObj = VideoWriter([movieName(1:end-4),'_2'],'mpeg-4');
% vidObj.FrameRate=movieFPS;
% vidObj.Quality=100;
% %     vidObj.LosslessCompression=1;
% open(vidObj);
% hh=waitbar(0,'please wait ...');
% cd(filePath)
% for ii=1:p
%                 waitbar(ii/p,hh, ...
%                  [num2str(100*ii/p,'%03.1f'),'% completed']);
%     file=fileName(ii).name;
%     img=differentTypeRead(file,fileType);
%     img2=newIntensity(img);
%     if get(handles.cb_eval,'value')
%         x=img2;
%         eval(get(handles.et_eval,'string')) ;
%         img2=x;
%     end
%     if size(img2,3)==3
% 
%         currFrame=im2frame(uint8(img2));
%                writeVideo(vidObj,currFrame);             
%     else
%     img3=uint8(img2(:,:));
%         currFrame=im2frame(img3,colorSelected);
%                writeVideo(vidObj,currFrame);        
%     end
% 
% end
% close(hh)
% close(vidObj);
% guidata(hObject,handles)
function et_fps_Callback(hObject, eventdata, handles)
% hObject    handle to et_fps (see GCBO)
% et_fps as text
%        str2double(get(hObject,'String')) returns contents of et_fps as a double
function et_fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_fps (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --------------------------------------------------------------------
function mp_pixelQuantity_Callback(hObject, eventdata, handles)
% hObject    handle to mp_pixelQuantity (see GCBO)
vv=0;
pointPst=handles.pointPst;
file=handles.file;
 img=differentTypeRead(file,handles.fileType);
disp(['%%%%%%%%%%%%%%%%%',datestr(now)])
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        bw=roipoly(img,cx,cy);
        vv=vv+1;
        if vv>1
            hold on;
        end
        disp(['group',num2str(ii),'ROI',num2str(jj),'pixel quantity:', num2str(sum(bw(:)))]);
    end
end
% --- Executes on button press in cb_unit16.
function cb_unit16_Callback(hObject, eventdata, handles)
function pushbutton99_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
cd(filePath)
img=differentTypeRead(file,fileType);
img=img(:,:,1);
if get(handles.cb_vertical,'value')==0
    img=img.';
end
[mm,nn]=size(img(:,:,1));
eval(get(handles.et_verticalLines,'string'));
p=length(yr);
Y=zeros(mm,p);
X=Y;
fig1=31;
fig2=32;
figure(fig1);close(fig1);figure(fig1);hold on;
figure(fig2);close(fig2);figure(fig2);hold on;
plotGroup={'r','g','b','r-*','g-*','b-*','r--','g--','b--'};
titleName='';
newImg=getCurrentImg(handles);
figure(fig2);imshow(uint8(newImg));hold on;
for ii=1:p
    yRange=yr{ii};
    y=mean(img(:,yRange),2);
    y=IOS_time_gui_filter(y,handles,'vertical');
    x=[1:mm].';
    figure(fig1);plot(x,y,plotGroup{ii});grid on;
    Y(:,ii)=y;
    X(:,ii)=x;
    titleName=[titleName,num2str(yRange(1)),' to ',num2str(yRange(end)),';'];
    figure(fig2);plot(0+1*y,x,plotGroup{ii})
end
title(titleName)
handles.vertical.Y=Y;
handles.vertical.X=X;
handles.vertical.Fig1=fig1;
handles.vertical.Fig2=fig2;
    guidata(hObject,handles)
% figure(31);plot(1:mm,y,1:mm,y,'ro'); grid on;
function et_verticalLines_Callback(hObject, eventdata, handles)
function et_verticalLines_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function result=filterName(result,handles)
handles=filterInfoRead2(handles);
filterMethod=handles.filterMethod;
sizNo=filterMethod.sizNo;
sigmaNo=filterMethod.sigmaNo;
nameNo=filterMethod.name;
if nameNo==1
    result=[result,'Gaus',num2str(sizNo),'_',num2str(sigmaNo)];
elseif nameNo==2
    result=[result,'Med',num2str(sizNo)];
elseif nameNo==3
    result=[result,'Wien',num2str(sizNo)];
elseif nameNo==0
end
function pushbutton100_Callback(hObject, eventdata, handles)
Y=handles.vertical.Y;
X=handles.vertical.X;
x=X(:,1);
eval(get(handles.et_verticalLines,'string'));
p=length(yr);

% y=handles.vertical.Y;
% x=handles.vertical.X;
fig1=handles.vertical.Fig1;
fig2=handles.vertical.Fig2;
dataSave=double([x,Y]);
cd(handles.filePath)
result='vertical';
result=filterName(result,handles);
resultFolder=fullfile(handles.filePath,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)

file=handles.file;
fileType=handles.fileType;

titleName='';

for ii=1:p
    yRange=yr{ii};
    titleName=[titleName,num2str(yRange(1),'%03d'),' to ',num2str(yRange(end),'%03d'),'_'];

end

if get(handles.cb_vertical,'value')
filenameTmp=[file(1:end-length(fileType)-1),'_Y',titleName,'.txt'];
else
filenameTmp=[file(1:end-length(fileType)-1),'_X',titleName,'.txt'];
end
[file,filePath]=uiputfile(filenameTmp);
cd(filePath)
figure(fig1)
saveas(gcf,file(1:end-4),'tiff')
figure(fig2)
saveas(gcf,[file(1:end-4),'com'],'tiff')
save(file,'dataSave','-ASCII')
cd(handles.filePath)
function pushbutton102_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
            cx=pointPst{1}{1}(:,1);
            cy=pointPst{1}{1}(:,2);
            cx(1:2)=1;
            %% deal with edge
            cx(cx==min(cx))=1;
            cx(cx==max(cx))=n;
bw=roipoly(img(:,:,1),cx,cy);
if get(handles.cb_sc,'value')
    yyPositionTmp=(1:m).';
    yyPosition=yyPositionTmp*ones(1,n);
    xxPositionTmp=(1:n).';
    xxPosition=xxPositionTmp*ones(1,m);
    xxPosition=xxPosition.';
    bw(:,1:8)=bw(:,9)*ones(1,8);
    bw(:,end-7:end)=bw(:,end-8)*ones(1,8);
    dist=sum(bw.*yyPosition)./sum(bw);
    dist=dist-dist(1);
    dist=IOS_time_gui_filter(dist,handles,'vertical');
    shift=dist;
     subShift=shift;
    subShift=subShift-floor(subShift);
    shift=floor(shift);
    shift=-shift;
%     dist=-dist;
    displacement=round(shift);
else
    %% binning size
    xBinning=6;
    yBinning=1;
    %% ROI
    xStart=1;
    xEnd=n;
    yStart=1;
    yEnd=m;
    %% verify whether ROI fits
    if mod(yEnd-yStart+1,yBinning)~=0
        yEnd=floor((yEnd-yStart+1)/yBinning)*yBinning-1+yStart;
    end
    if mod(xEnd-xStart+1,xBinning)~=0
        xEnd=floor((xEnd-xStart+1)/xBinning)*xBinning-1+xStart;
    end
    new_xsize=(xEnd-xStart+1)/xBinning; % image size after pixel binning
    new_ysize=(yEnd-yStart+1)/yBinning; % image size after pixel binning
    imgBin=image_binning(img,xBinning,yBinning,xStart,yStart,new_xsize,new_ysize);
    bwBin=image_binning(bw,xBinning,yBinning,xStart,yStart,new_xsize,new_ysize);
    bwBin=bwBin>0;
    %% find displacement
    shiftBin=zeros(new_xsize,1);
    bwBin=double(bwBin);
    dataTmp=imgBin.*bwBin;
    for ii=1:new_xsize
        data=dataTmp(:,ii);
        maxData=max(data);
         if maxData~=0
             I=find(data==maxData);
             shiftBin(ii)=I(1);
         end
    end
    shiftBin(1:3)=shiftBin(4);
    shiftBin=IOS_time_gui_filter(shiftBin,handles,'vertical');
    shiftBin=round(shiftBin);
    displacement=zeros(n,1);
    for ii=1:new_xsize-1
        displacement((ii-1)*xBinning+1:ii*xBinning)=shiftBin(ii);
    end
    ii=new_xsize;
    displacement((ii-1)*xBinning+1:end)=shiftBin(ii);
    displacement=-displacement;
end
%%
figure(21); plot(1:n,displacement)
minY=min(displacement)-1;
maxY=max(displacement)+1;
if minY>0
    minY=-1;
end
if strcmp(fileType,'mat')==0
imgInfo=imfinfo(fileName(1).name);
uint8Flag=imgInfo.BitDepth;
end
imgNew=zeros(m-minY+maxY,n,['uint',num2str(uint8Flag)]);
for ii=1:n
    imgNew(-minY+displacement(ii):-minY+displacement(ii)+m-1,ii)=img(:,ii);
end
if get(handles.cb_sc,'value')
    imgNew=img_ROI_registerSubpixel(imgNew,subShift,'bilinear');
    imgNew=eval(['uint',num2str(uint8Flag),'(imgNew)']);
end
result='deCurve';
resultPath=fullfile(handles.filePath,result);
cd(handles.filePath)
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
[fileTmp,resultPath]=uiputfile(file);
cd(resultPath)
imwrite(imgNew,fileTmp)
function imgROI2=img_ROI_registerSubpixel(img,subShift,method)
[m,n]=size(img);
m=m-1;
n=n-1;
imgROI2=zeros(m,n);
img=double(img);
deltaY=subShift(1);
deltaX=subShift(2);
for ii=1:m
    for jj=1:n
        if strcmp(method,'bilinear') || strcmp(method,'Bilinear')
        %% bilinear interpolation
        a=img(ii,jj); b=img(ii,jj+1);
        c=img(ii+1,jj+1);d=img(ii+1,jj);
        b1=a;
        b2=b-a;
        b3=d-a;
        b4=a+c-b-d;
        imgROI2(ii,jj)=b1+b2*deltaX+b3*deltaY+b4*deltaX*deltaY;
        end
        if strcmp(method,'linear') || strcmp(method,'Linear')
            imgROI2(ii,jj)=[1-deltaX,deltaX]*(img(ii:ii+1,jj:jj+1)).'*[1-deltaY,deltaY].';
        end
    end
end
function pushbutton103_Callback(hObject, eventdata, handles)
if get(handles.cb_rotation,'value')
    cd(handles.filePath)
    fileName=handles.fileName;
    p=length(fileName);
    cd(handles.filePath)
    result='rotate';
    resultFolder=fullfile(handles.filePath,result);
    if exist(resultFolder)==7
    else
        mkdir(result)
    end
    cd(resultFolder)
    [file,resultPath]=uiputfile(handles.file);
    h_wait=waitbar(0,'please wait');
    if strcmp(handles.fileType,'mat')~=1
        for ii=1:p
            waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
            file=fileName(ii).name;
            img=imread(fullfile(handles.filePath,file));
            eval(get(handles.et_rotateAngle,'string'));
            %img= imrotate(img,degrees,'nearest','crop');
            imwrite(imgNew,fullfile(resultPath,file),handles.fileType);
        end
    else
        for ii=1:p
            waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
            file=fileName(ii).name;
            cd(handles.filePath)
            img=differentTypeRead(file,handles.fileType);
            eval(get(handles.et_rotateAngle,'string'));
            %img= imrotate(img,degrees,'nearest','crop');
                 cd(resultPath)
            save(file,'imgNew');
        end
    end
    close(h_wait)
end
function et_rotateAngle_Callback(hObject, eventdata, handles)
function et_rotateAngle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton100_CreateFcn(hObject, eventdata, handles)
% --------------------------------------------------------------------
function me_background_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function me_G1ROI1_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
cx=pointPst{1}{1}(:,1);
cy=pointPst{1}{1}(:,2);
x1=min(cx(:));
x2=max(cx(:));
y1=min(cy(:));
y2=max(cy(:));
x1=round(x1);x2=round(x2);y1=round(y1);y2=round(y2);
imgROI=img(y1:y2,x1:x2,1);
background=zeros(m,1);
background(y1:y2)=mean(imgROI,2);
background=IOS_time_gui_filter(background,handles,'vertical');
figure(33);plot(1:length(background),background);title('background');xlabel('y')
handles.background=background*ones(1,n);
guidata(hObject,handles);
 %   uiwait(msgbox('DONE'));
% --------------------------------------------------------------------
function me_backgroundRun_Callback(hObject, eventdata, handles)
background=handles.background;
if strcmp(handles.fileType,'mat')~=1
    background=eval([handles.imgDepth,'(handles.background)']);
end
    fileName=handles.fileName;
    p=length(fileName);
    [file,resultPath]=uiputfile(handles.file);
    h_wait=waitbar(0,'please wait');
    if strcmp(handles.fileType,'mat')~=1
        for ii=1:p
            waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
            file=fileName(ii).name;
            img=imread(fullfile(handles.filePath,file));
            %img= imrotate(img,degrees,'nearest','crop');
            imwrite(img-background,fullfile(resultPath,file),handles.fileType,'compression','none');
        end
    else
        for ii=1:p
            waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
            file=fileName(ii).name;
            cd(handles.filePath)
            img=differentTypeRead(file,handles.fileType);
            imgNew=img-background;
            %img= imrotate(img,degrees,'nearest','crop');
                 cd(resultPath)
            save(file,'imgNew');
        end
    end
    close(h_wait)
function et_db_y0_Callback(hObject, eventdata, handles)
function et_db_y0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_db_dy_Callback(hObject, eventdata, handles)
function et_db_dy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function db_I=db_getI(a,b)
a_b=abs(a-b);
[Y,db_I]=min(a_b);
% --- Executes on button press in et_db_single.
function et_db_single_Callback(hObject, eventdata, handles)
% hObject    handle to et_db_single (see GCBO)
y0=eval(get(handles.et_db_y0,'string'));
dy=eval(get(handles.et_db_dy,'string'));
file=handles.file;
handles=filterInfoRead2_spatial(handles);%img=imgFilter(img,filterMethod2)
filterMethod2=handles.filterMethod2;
cd(handles.filePath);
img=differentTypeRead(file,handles.fileType);
[mm,nn,n_trash]=size(img);
if get(handles.pp_db_saw,'value')==1
    if n_trash==1
        a=img(y0-dy:y0,:);b=img(y0:y0+dy,:);
        db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew=[img(1:y0-dy-1,:);c;img(y0+dy+1:end,:)];
        imgNew=imgFilter(imgNew,filterMethod2);
    elseif n_trash==3
        %r
         a=img(y0-dy:y0,:,1);b=img(y0:y0+dy,:,1);
        db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew1=[img(1:y0-dy-1,:,1);c;img(y0+dy+1:end,:,1)];
         %g
         a=img(y0-dy:y0,:,2);b=img(y0:y0+dy,:,2);
    %     db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew2=[img(1:y0-dy-1,:,2);c;img(y0+dy+1:end,:,2)];
             %b
         a=img(y0-dy:y0,:,3);b=img(y0:y0+dy,:,3);
    %     db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew3=[img(1:y0-dy-1,:,3);c;img(y0+dy+1:end,:,3)];
        imgNew=zeros(size(imgNew3,1),size(imgNew2,2),3);
        imgNew1=imgFilter(imgNew1,filterMethod2);
        imgNew2=imgFilter(imgNew2,filterMethod2);
        imgNew3=imgFilter(imgNew3,filterMethod2);
        imgNew(:,:,1)=imgNew1;
        imgNew(:,:,2)=imgNew2;
        imgNew(:,:,3)=imgNew3;
    end
elseif get(handles.pp_db_saw,'value')==2
        a=img(y0-dy:y0+dy,:,:);
        db_ratio=eval(get(handles.et_db_ratio,'string'));
        numrows=round(size(a,1)*db_ratio);
        numcols=size(a,2);
        c=imresize(a,[numrows,numcols]);
        if n_trash==1
             imgNew=[img(1:y0-dy-1,:,1);c(:,:,1);img(y0+dy+1:end,:,1)];
        elseif n_trash==3
            imgNew1=[img(1:y0-dy-1,:,1);c(:,:,1);img(y0+dy+1:end,:,1)];
            imgNew2=[img(1:y0-dy-1,:,2);c(:,:,2);img(y0+dy+1:end,:,2)];
            imgNew3=[img(1:y0-dy-1,:,3);c(:,:,3);img(y0+dy+1:end,:,3)];
            imgNew=zeros(size(imgNew3,1),size(imgNew2,2),3);
         imgNew1=imgFilter(imgNew1,filterMethod2);
        imgNew2=imgFilter(imgNew2,filterMethod2);
        imgNew3=imgFilter(imgNew3,filterMethod2);
            imgNew(:,:,1)=imgNew1;
            imgNew(:,:,2)=imgNew2;
            imgNew(:,:,3)=imgNew3;
        end
end
    %%
    contents = cellstr(get(handles.pm_colormap_raw,'String'));
    colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
    colorSelected=eval(colorSelectedTmp{1});
    handles.newIntensity=inline(get(handles.et_inline,'string'));
    newImg=handles.newIntensity(imgNew);
    d3=size(newImg,3);
        %set(h_img,'CData',uint16(img))
        if d3==3
            figure(31);imshow(uint8(imgNew))
        elseif d3==1
            newImgColor=ind2rgb(uint8(newImg),colorSelected);
           figure(31);imshow(uint8(newImgColor))
          %  set(h_img,'CData',uint8(newImg))
        end
function et_db_runAll_Callback(hObject, eventdata, handles)
y0=eval(get(handles.et_db_y0,'string'));y0=round(y0);
dy=eval(get(handles.et_db_dy,'string'));dy=round(dy);
file=handles.file;
handles=filterInfoRead2_spatial(handles);%img=imgFilter(img,filterMethod2)
filterMethod2=handles.filterMethod2;
fileName=handles.fileName;
p=length(fileName);
cd(handles.filePath)
[result,resultPath]=uiputfile(file);
for ii=1:p
    file=handles.fileName(ii).name;
    cd(handles.filePath);
    img=differentTypeRead(file,handles.fileType);
    [mm,nn,n_trash]=size(img);
if get(handles.pp_db_saw,'value')==1
    if n_trash==1
        a=img(y0-dy:y0,:);b=img(y0:y0+dy,:);
        db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew=[img(1:y0-dy-1,:);c;img(y0+dy+1:end,:)];
        imgNew=imgFilter(imgNew,filterMethod2);
    elseif n_trash==3
        %r
         a=img(y0-dy:y0,:,1);b=img(y0:y0+dy,:,1);
        db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew1=[img(1:y0-dy-1,:,1);c;img(y0+dy+1:end,:,1)];
         %g
         a=img(y0-dy:y0,:,2);b=img(y0:y0+dy,:,2);
    %     db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew2=[img(1:y0-dy-1,:,2);c;img(y0+dy+1:end,:,2)];
             %b
         a=img(y0-dy:y0,:,3);b=img(y0:y0+dy,:,3);
    %     db_I=db_getI(a,b);
        c=a;
        for jj=1:nn
            c(db_I(jj)+1:end,jj)=b(db_I(jj)+1:end,jj);
        end
        imgNew3=[img(1:y0-dy-1,:,3);c;img(y0+dy+1:end,:,3)];
        imgNew=zeros(size(imgNew3,1),size(imgNew2,2),3);
        imgNew1=imgFilter(imgNew1,filterMethod2);
        imgNew2=imgFilter(imgNew2,filterMethod2);
        imgNew3=imgFilter(imgNew3,filterMethod2);
        imgNew(:,:,1)=imgNew1;
        imgNew(:,:,2)=imgNew2;
        imgNew(:,:,3)=imgNew3;
    end
elseif get(handles.pp_db_saw,'value')==2
        a=img(y0-dy:y0+dy,:,:);
        db_ratio=eval(get(handles.et_db_ratio,'string'));
        numrows=round(size(a,1)*db_ratio);
        numcols=size(a,2);
        c=imresize(a,[numrows,numcols]);
        if n_trash==1
             imgNew=[img(1:y0-dy-1,:,1);c(:,:,1);img(y0+dy+1:end,:,1)];
        elseif n_trash==3
            imgNew1=[img(1:y0-dy-1,:,1);c(:,:,1);img(y0+dy+1:end,:,1)];
            imgNew2=[img(1:y0-dy-1,:,2);c(:,:,2);img(y0+dy+1:end,:,2)];
            imgNew3=[img(1:y0-dy-1,:,3);c(:,:,3);img(y0+dy+1:end,:,3)];
            imgNew=zeros(size(imgNew3,1),size(imgNew2,2),3);
         imgNew1=imgFilter(imgNew1,filterMethod2);
        imgNew2=imgFilter(imgNew2,filterMethod2);
        imgNew3=imgFilter(imgNew3,filterMethod2);
            imgNew(:,:,1)=imgNew1;
            imgNew(:,:,2)=imgNew2;
            imgNew(:,:,3)=imgNew3;
        end
end
cd(resultPath)
if strcmp(handles.fileType,'mat')~=1
    if n_trash==3
        imwrite(uint8(imgNew),handles.fileName(ii).name);
    else
        imwrite(eval([handles.imgDepth,'(imgNew)']),handles.fileName(ii).name);
    end
else
    save(handles.fileName(ii).name,'imgNew')
end
end
% --- Executes on selection change in pp_db_saw.
function pp_db_saw_Callback(hObject, eventdata, handles)
function pp_db_saw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_db_ratio_Callback(hObject, eventdata, handles)
function et_db_ratio_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function img=imgFilter(img,filterMethod2)
sizNo=filterMethod2.sizNo;
sigmaNo=filterMethod2.sigmaNo;
nameNo=filterMethod2.name;
%% filterImg
if nameNo==0
    %imgAvg=imgAvg;
elseif nameNo==1
    w=fspecial('gaussian',[sizNo,sizNo],sigmaNo);
    img=imfilter(img(:,:,1),w,'replicate');
elseif nameNo==2
    img=medfilt2(img(:,:,1),[sizNo,sizNo]);
elseif nameNo==3
        img=wiener2(img(:,:,1),[sizNo,sizNo]);
elseif nameNo==4
    imgAvg2=edge(img(:,:,1),'canny',[.05,.4],sigmaNo);
    img=single(imgAvg2);
end
% --------------------------------------------------------------------
function me_backgroundRun_divide_Callback(hObject, eventdata, handles)
cd(handles.filePath)
background=handles.background;
background(background<=1)=1;
background=background/max(background(:,1));
fileName=handles.fileName;
fileType=handles.fileType;
p=length(fileName);
         saveName=[handles.file(1:end-length(fileType)-1),'.mat'];
[file,resultPath]=uiputfile(saveName);
h_wait=waitbar(0,'please wait');
for ii=1:p
    waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
    file=fileName(ii).name;
    cd(handles.filePath)
    img=differentTypeRead(file,handles.fileType);
%     img2=img;
    imgNew=img./background;
%     temp=img2>max(0.4*img2(:));
%     imgNew(temp)=img2(temp);
    %img= imrotate(img,degrees,'nearest','crop');
         cd(resultPath)
         saveName=[file(1:end-length(fileType)-1),'.mat'];
    save(saveName,'imgNew');
end
    close(h_wait)
% --- Executes on button press in cb_saveMap.
function cb_saveMap_Callback(hObject, eventdata, handles)
function pushbutton108_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton108 (see GCBO)
[file,filePath]=uigetfile('*.*','select the map');
dotNO=strfind(file,'.');dotNO=dotNO(end);
fileType=file(dotNO+1:end);
cd(filePath)
img=differentTypeRead(file,fileType);
handles.vesselMask=img;
guidata(hObject,handles)
% --- Executes on button press in cb_vessel_mouse.
function cb_vessel_mouse_Callback(hObject, eventdata, handles)
function et_vessel_bg_Callback(hObject, eventdata, handles)
function et_vessel_bg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton109_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton109 (see GCBO)
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
aa=handles.vesselMask;
vesselMask=aa>max(aa(:))*.2;
vesselMask=logical(aa);
vesselMask=~vesselMask;
cd(filePath)
img=differentTypeRead(file,fileType);
if get(handles.cb_vessel_mouse,'value')
    contents = cellstr(get(handles.pm_colormap_raw,'String'));
    colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
    colorSelected=eval(colorSelectedTmp{1});
%     cd(handles.filePath)
    handles.newIntensity=inline(get(handles.et_inline,'string'));
    newImg=handles.newIntensity(img);
    d3=size(newImg,3);
    figure(39);
         if d3==3
            imshow(uint8(newImg))
        elseif d3==1
            newImgColor=ind2rgb(uint8(newImg),colorSelected);
            imshow(newImgColor)
           %imshow(uint8(newImg),colorSelected)
         end
        [c,r,p]=impixel;
        background=img(c(end),r(end),:);
else
    background=eval(get(handles.et_vessel_bg,'string'));
end
[mm,nn,nn_trash]=size(img);
cd(filePath)
[file2,resultPath]=uiputfile(file,'save files');
fileName=handles.fileName;p=length(fileName);
h_wait=waitbar(0,'please wait');
for ii=1:p
    waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
    file=fileName(ii).name;
    cd(handles.filePath)
    img=differentTypeRead(file,handles.fileType);
    if nn_trash==3
        imgtmp=img(:,:,1);imgtmp(vesselMask)=background(1);img(:,:,1)=imgtmp;
        imgtmp=img(:,:,2);imgtmp(vesselMask)=background(2);img(:,:,2)=imgtmp;
        imgtmp=img(:,:,3);imgtmp(vesselMask)=background(3);img(:,:,3)=imgtmp;
    else
        img(vesselMask)=background(1);
    end
    cd(resultPath)
%     img2=img;
%     temp=img2>max(0.4*img2(:));
%     imgNew(temp)=img2(temp);
    %img= imrotate(img,degrees,'nearest','crop');
    if strcmp(fileType,'mat')
         saveName=[file(1:end-length(fileType)-1),'.mat'];
    save(saveName,'img');
    else
        imwrite(eval([handles.imgDepth,'(img)']),file)
    end
end
close(h_wait)
function pushbutton110_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton110 (see GCBO)
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
cd(filePath)
img0=differentTypeRead(file,fileType);
[m,n,n_trash]=size(img0);
fun=get(handles.et_imresize,'string');
result='resize';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(result);
end
cd(result)
[result,resultPath]=uiputfile(file);
imgDepth=handles.imgDepth;
h_wait=waitbar(0,'please wait');
for ii=1:p
    waitbar(ii/p,h_wait,num2str(100*ii/p,'%03.1f'));
    file=fileName(ii).name;
    cd(filePath)
    a_tmp=differentTypeRead(file,fileType);
    if n_trash==1
        a=a_tmp;
        eval(fun);
        img=b;
    elseif n_trash==3
        a=a_tmp(:,:,1);
        eval(fun);
        img(:,:,1)=b;
        a=a_tmp(:,:,2);
        eval(fun);
        img(:,:,2)=b;
         a=a_tmp(:,:,3);
        eval(fun);
        img(:,:,3)=b;
    end
    cd(resultPath)
    differentTypeWrite(img,file, fileType,imgDepth);
end
close(h_wait)
function  differentTypeWrite(img,file, fileType,imgDepth)
if strcmp(fileType,'mat')
    save(file,'img');
else
     imwrite(eval([imgDepth,'(img)']),file)
end
function et_imresize_Callback(hObject, eventdata, handles)
function et_imresize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --------------------------------------------------------------------
function me_comb_Callback(hObject, eventdata, handles)
function pushbutton111_Callback(hObject, eventdata, handles)
rg=eval(get(handles.et_gr,'string'));
handles.pointPst{rg(1)}(rg(2))=[];
delete(handles.pointPstH{rg(1)}(rg(2)))
handles.pointPstH{rg(1)}(rg(2))=[];
handles=GroupROITag(handles);
guidata(hObject,handles)
function et_gr_Callback(hObject, eventdata, handles)
function et_gr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --------------------------------------------------------------------
function me_half_half_Callback(hObject, eventdata, handles)
cd(handles.filePath)
%% imread 1
file1=handles.file;
filePath1=handles.filePath;
% [file1,filePath1]=uigetfile('*.*','select 1st series of images');
% cd(filePath1);
indexNO=strfind(file1,'.');
fileType1=file1(indexNO(end)+1:end);
fileName1=dir(['*.',fileType1]);
p=length(fileName1);
img=imread(file1);
[m1,n1,q1]=size(img);
img3=zeros(size(img));
imgInfo=imfinfo(fileName1(1).name);
uint8Flag1=imgInfo.BitDepth;
if uint8Flag1==24
    uint8Flag1=8;
elseif uint8Flag1==16*3
    uint8Flag1=16;
end
%% imread 2
[file2,filePath2]=uigetfile('*.*','select 2n series of images');
cd(filePath2);
indexNO=strfind(file2,'.');
fileType2=file2(indexNO(end)+1:end);
fileName2=dir(['*.',fileType2]);
img=imread(file2);
[m2,n2,q2]=size(img);
cd(filePath1)
result='combined';
resultFolder=fullfile(filePath1,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
[result,resultPath]=uiputfile(file1,'choose a folder to save');
%%
 prompt={'ratio:'};
   name='Input for Peaks function';
   numlines=1;
   defaultanswer={'[0.5,0.5]'};
    options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
   answer=inputdlg(prompt,name,numlines,defaultanswer,options);
   ratio=eval(answer{1});
 %  ratio=ratio/sum(ratio);
%%
h_wait=waitbar(0,'wait');
if q1==3 && q2==1
    for ii=1:p
            waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
        cd(filePath1)
        file1=fileName1(ii).name;
        img1=double(imread(file1));
        cd(filePath2)
        file2=[file1(1:end-length(fileType1)),fileType2];
        img2=double(imread(file2));
        img3(:,:,1)=ratio(1)*img1(:,:,1)+ratio(2)*img2;
        img3(:,:,2)=ratio(1)*img1(:,:,2)+ratio(2)*img2;
        img3(:,:,3)=ratio(1)*img1(:,:,3)+ratio(2)*img2;
        cd(resultPath)
        differentTypeWrite(img3,file1, handles.fileType,handles.imgDepth);
%         imwrite(img1,file1);
    end
elseif q1==1 && q2==1
     for ii=1:p
            waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
        cd(filePath1)
        file1=fileName1(ii).name;
        img1=differentTypeRead(file1,fileType1);
        cd(filePath2)
        file2=fileName2(ii).name;
        img2=differentTypeRead(file2,fileType2);
        img3(:,:,1)=ratio(1)*img1(:,:,1)+ratio(2)*img2;
        cd(resultPath)
        differentTypeWrite(img3(:,:,1),file1, handles.fileType,handles.imgDepth);
%         imwrite(img1,file1);
    end
end
close(h_wait)
% --------------------------------------------------------------------
function me_partialCombine_Callback(hObject, eventdata, handles)
% hObject    handle to me_partialCombine (see GCBO)
 uiwait(msgbox('make sure the current folder is the one with color images','Title','sure'));
   [file1,filePath1]=uigetfile('*.*','select map that indicates the pixels to be left');
    cd(filePath1)
file=handles.file;
filePath=handles.filePath;
  %  [file,filePath]=uigetfile('*.*','select IOS image files');
    [file2,filePath2]=uigetfile('*.*','select RAW image files');
    result='superimpose';
%% imread files
cd(filePath1);
indexNO1=strfind(file1,'.');
fileType1=file1(indexNO1(end)+1:end);
aa=differentTypeRead(file1,fileType1);
% vesselMask=aa>max(aa(:))*.2;
vesselMask=logical(aa);
vesselMask=~vesselMask;
cd(filePath);
indexNO=strfind(file,'.');
fileType=file(indexNO(end)+1:end);
fileName=dir(['*.*',fileType]);
p=length(fileName);
indexNO=strfind(file2,'.');
fileType2=file2(indexNO(end)+1:end);
%% directory for results
resultPath=fullfile(filePath,result);
cd(filePath)
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
%%
h_wait=waitbar(0,'wait');
for ii=1:p
      waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
    cd(filePath)
    file=fileName(ii).name;
    imgIOS=imread(file);
    cd(filePath2)
    imgRaw=imread([file(1:end-length(fileType)),fileType2]);
    cd(resultPath)
    img1=imgIOS(:,:,1);
    img2=imgIOS(:,:,2);
    img3=imgIOS(:,:,3);
    img1(vesselMask)=imgRaw(vesselMask);
    img2(vesselMask)=imgRaw(vesselMask);
    img3(vesselMask)=imgRaw(vesselMask);
    imgIOS(:,:,1)=img1;
    imgIOS(:,:,2)=img2;
    imgIOS(:,:,3)=img3;
    imwrite(uint8(imgIOS),file)
end
close(h_wait)
function et_shiftEst_fn_Callback(hObject, eventdata, handles)
function et_shiftEst_fn_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_shiftEst.
function pb_shiftEst_Callback(hObject, eventdata, handles)

pointPst=handles.pointPst;

p=length(pointPst);
cx=cell(p,1);
cy=cell(p,1);
%%
% pointPst{1}(2)=[]
%%
for ii=1:p
    a=pointPst{ii};
    q=length(a);
    cx{ii}=cell(q,1);
    cy{ii}=cell(q,1);
    for jj=1:q
        cx{ii}{jj}=a{jj}(:,1);
        cy{ii}{jj}=a{jj}(:,2);
        
        
    end
    
end
file=handles.file;
filePath=handles.filePath;
% [file,filePath]=uigetfile('*.*','select raw');
folderName=get(handles.et_shiftEst_fn,'string');
wt=eval(get(handles.et_shiftEst_wt,'string'));
scaleWeight=wt(1);
thresh=wt(2);
othersH.bigger=wt(3);

handles=filterInfoRead2(handles);
filterMethod=handles.filterMethod;
sizNo=filterMethod.sizNo;
sigmaNo=filterMethod.sigmaNo;
nameNo=filterMethod.name;
othersH.filterMethod=filterMethod;

for ii=1:p
    if ~isempty(pointPst{ii})
        spotsNumber=length(cx{ii});
        [file,filePath,shift,y_saver2,x_saver2]= ...
            image_shiftEstimate_centroid_gui(file,filePath,handles.fileName, ...
            spotsNumber,cy{ii},cx{ii},['group',num2str(ii)],folderName, ...
            scaleWeight,thresh,othersH);        
    end

end
result0=folderName;
resultPath0=fullfile(filePath,result0);
cd(resultPath0)
if exist(resultPath0)==7
else
    mkdir(result0)
end
pointPst=handles.pointPst;
colorTmp=handles.colorTmp;
save(['group',num2str(ii),'parameterMat.mat'],'pointPst','colorTmp')
handles.shift=shift;

guidata(hObject,handles)
% pointPst=handles.pointPst;
% p=length(pointPst);
% cx=cell(p,1);
% cy=cell(p,1);
% %%
% % pointPst{1}(2)=[]
% %%
% for ii=1:p
%     a=pointPst{ii};
%     q=length(a);
%     cx{ii}=cell(q,1);
%     cy{ii}=cell(q,1);
%     for jj=1:q
%         cx{ii}{jj}=a{jj}(:,1);
%         cy{ii}{jj}=a{jj}(:,2);
%     end
% end
% file=handles.file;
% filePath=handles.filePath;
% % [file,filePath]=uigetfile('*.*','select raw');
% folderName=get(handles.et_shiftEst_fn,'string');
% wt=eval(get(handles.et_shiftEst_wt,'string'));
% scaleWeight=wt(1);
% thresh=wt(2);
% othersH.bigger=wt(3);
% for ii=1:p
%     spotsNumber=length(cx{ii});
%     [file,filePath,shift,y_saver2,x_saver2]= ...
%         image_shiftEstimate_centroid_gui(file,filePath, ...
%         spotsNumber,cy{ii},cx{ii},['group',num2str(ii)],folderName, ...
%         scaleWeight,thresh,othersH)
% end
function [file,filePath,shift,y_saver2,x_saver2]=image_shiftEstimate_centroid_gui(file,filePath,fileName,spotsNumber,y_saver2,x_saver2,folderGroup,folder0,scaleWeight,thresh, othersH)
aa=[1,2,3,4,11, 5];
for ii=1:length(aa)
    figure(aa(ii))
    close(aa(ii))
end
if nargin<2
    [file,filePath]=uigetfile('*.*');
end
if nargin<3
    spotsNumber=1;
end
dotNO=strfind(file,'.');
fileType=file(dotNO(end)+1:end);
cd(filePath)
img=differentTypeRead(file,fileType);
figure(1); h_axis=gca;
% imgHandle=imfinfo(file);
% BitDepth=imgHandle.BitDepth;
% imageShow(h_axis,BitDepth,img);
imshow(im2uint8(mat2gray(double(img))))
hold on;
ROIii={'r-','g-','b-','r:','g:','b:','r.','g.','b.','r-.','g-.','b-.','*-','*--','*:','*-.'};
displacement=cell(spotsNumber,1);
margin=0;
%% define ROIs
if nargin<5
%     x_saver=zeros(spotsNumber,2);
%     y_saver=zeros(spotsNumber,2);
     x_saver2=cell(spotsNumber,1);
     y_saver2=cell(spotsNumber,1);
    for ii=1:spotsNumber
         uiwait(msgbox(['select ROI ',num2str(ii)]));
         axes(h_axis);
         [cx,cy,cc,x_saverTmp2,y_saverTmp2]=improfile;
         x_saverTmp2=round(x_saverTmp2);
         y_saverTmp2=round(y_saverTmp2);
         c1=x_saverTmp2;r1=y_saverTmp2;
%          [c1,r1,p_tmp]=impixel;
         c2=sort(c1);
         r2=sort(r1);
        x1=c2(1)-margin;x2=c2(end)+margin;
        y1=r2(1)-margin;y2=r2(end)+margin;
%          x1=118;x2=131;y1=50;y2=61;
%         x_saver(ii,:)=[x1,x2];
%         y_saver(ii,:)=[y1,y2];
         plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],ROIii{ii})
         hold on;
         x_temp=[x_saverTmp2.',x_saverTmp2(1)];
         y_temp=[y_saverTmp2.',y_saverTmp2(1)];
         plot(x_temp,y_temp,'r*-')
         x_saver2{ii}=x_saverTmp2;
         y_saver2{ii}=y_saverTmp2;
    end
end
x_saver=zeros(spotsNumber,2);
y_saver=zeros(spotsNumber,2);
for ii=1:spotsNumber
    x_saver2{ii}=round(x_saver2{ii});
    y_saver2{ii}=round(y_saver2{ii});
    [cx,cy,cc,x_saverTmp2,y_saverTmp2]=improfile(img,x_saver2{ii},y_saver2{ii});
    c1=x_saverTmp2;r1=y_saverTmp2;
    c2=sort(c1);
    r2=sort(r1);
    x1=c2(1)-margin;x2=c2(end)+margin;
    y1=r2(1)-margin;y2=r2(end)+margin;
    x_saver(ii,:)=[x1,x2];
    y_saver(ii,:)=[y1,y2];
         plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],ROIii{ii},'lineWidth',2)
         hold on;
         x_temp=[x_saverTmp2.',x_saverTmp2(1)];
         y_temp=[y_saverTmp2.',y_saverTmp2(1)];
         plot(x_temp,y_temp,'r-','lineWidth',2)
end
%%
filterMethod=othersH.filterMethod;
% filterMethod.name=2; % 0: none; 1: gaussian; 2: median; 3: Wiener
% filterMethod.sizNo=5;
% filterMethod.sigmaNo=3; % no used for median & wiener filters
%% create folder
cd(filePath)
if nargin<=5
    folderGroup=['group1'];
end
if nargin<=6
    folder0=['shift result'];
end
result=folderGroup;
result0=folder0;
resultPath0=fullfile(filePath,result0);
cd(filePath)
if exist(resultPath0)==7
else
    mkdir(result0)
end
resultPath=fullfile(resultPath0,result);
cd(resultPath0)
if exist(resultPath)==7
else
    mkdir(result)
end
%% find displacements
energyTmp=cell(spotsNumber,1);
areaTmp=cell(spotsNumber,1);
for ii=1:spotsNumber
    x1=x_saver(ii,1);x2=x_saver(ii,2);
    y1=y_saver(ii,1);y2=y_saver(ii,2);
    [centroids,energyTmp{ii},areaTmp{ii}]= ...
        centroidsFinder(file,filePath,fileName,x1,x2,y1,y2,x_saver2{ii},y_saver2{ii}, ...
        filterMethod,resultPath,ii,scaleWeight,thresh,othersH);
    p=size(centroids,1);
    relativeCentroids=centroids-ones(p,1)*centroids(1,:);
    displacement{ii}=relativeCentroids;
end
%% average displacements
shift=zeros(p,2);
energy=zeros(p,1);
area=zeros(p,1);
for ii=1:spotsNumber
    shift=shift+displacement{ii};
    energy=energy+energyTmp{ii};
    area=area+areaTmp{ii};
end
shift=shift/spotsNumber;
energy=energy/spotsNumber;
area=area/spotsNumber;
% shift=shift*2;
figure(2); h1=subplot(2,1,1); plot(1:p,shift(:,2),'r+-'); xlabel('avg horizontal shift')
h2=subplot(2,1,2); plot(1:p,shift(:,1),'r+-'); xlabel('avg vertical shift');
%% filter shifts
shift(:,2)=IOS_time_gui_filter_shift(shift(:,2),filterMethod);
axes(h1); hold on; plot(1:p,shift(:,2),'b+-'); legend({'raw','filtered'});
shift(:,1)=IOS_time_gui_filter_shift(shift(:,1),filterMethod);
axes(h2); hold on; plot(1:p,shift(:,1),'b+-'); legend({'raw','filtered'});
% img_register(file,filePath,shift);
energy=IOS_time_gui_filter_shift(energy,filterMethod);
area=IOS_time_gui_filter_shift(area,filterMethod);
figure(3);subplot(2,1,1);plot(1:p,energy,'r*-');title('energy')
subplot(2,1,2);plot(1:p,area,'r*-');title('area')
%%
cd(resultPath)
result='result_all';
result1=['shift_indvdl'];
resultPath1=fullfile(resultPath,result1);
cd(resultPath)
if exist(resultPath1)==7
else
    mkdir(result1)
end
cd(resultPath1)
displacement_raw=displacement;
for ii=1:spotsNumber
    displacement{ii}=IOS_time_gui_filter_shift(displacement{ii},filterMethod);
 areaTmp{ii}=IOS_time_gui_filter_shift(areaTmp{ii},filterMethod);
end
figure(5);
set(gcf,'position',[  1           1        1064         410])
h1=axes;
set(h1,'position',[ 0.1281    0.5814    0.7750    0.3412])
hold on;
h2=axes;
set(h2,'position',[0.1300    0.1100    0.7750    0.3412])
hold on;
% colorTmp=[1,0,0;
%           0,1,0;
%           0,0,1;
%           0.502,0,1;
%           0,1,1;
%           1,0,1;
%           0.3,0.3,0.7;
%           0.5,1/2,1/2];
legendName=cell(spotsNumber,1);
for ii=1:spotsNumber
    a=displacement{ii};
    b=displacement_raw{ii};
    p=length(a);
    t=[1:p].';
    shiftSave=[t,a,b];
    save(['cell ',num2str(ii,'%02d'),'_y_x_yraw_xraw.txt'],'shiftSave','-ASCII')
    axes(h1) %#ok<*LAXES>
    plot(t,a(:,2),ROIii{ii})
    axes(h2)
    plot(t,a(:,1),ROIii{ii})
    legendName{ii}=['spots',num2str(ii)];
end
cd ..
axes(h1);legend(legendName);xlabel('horizontal+:right');
axes(h2);legend(legendName);ylabel('vertical(+:down)')
saveas(gcf,'shift.tiff')
% save([result,'_area.txt'],'area','-ASCII')
% save([result,'_energy.txt'],'energy','-ASCII')
save([result,'_Avgshift.txt'],'shift','-ASCII')
save(result,'x_saver2','y_saver2','shift','energyTmp','energy','areaTmp','area','displacement')
figure(2);
saveas(gcf,'avg.tiff')
figure(1);
saveas(gcf,'cellsMap.tiff')
function [centroids,varargout]=centroidsFinder(file,filePath,fileName,x1,x2,y1,y2,x_saver2,y_saver2,filterMethod,resultPath0,cellNO,scaleWeight,threshold,othersH)
%% image imread
% [file,filePath]=uigetfile('*.*');
dotNO=strfind(file,'.');
fileType=file(dotNO(end)+1:end);
cd(filePath)
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
img=differentTypeRead(file,fileType);
%% ROI choose
% figure;
 imgHandle=imfinfo(file);
 BitDepth=imgHandle.BitDepth;
% imageShow(gca,BitDepth,img);
%
% % [c1,r1,p_tmp]=impixel;
% % c2=sort(c1);
% % r2=sort(r1);
% %
% % x1=c2(1);x2=c2(end);
% % y1=r2(1);y2=r2(end);
% x1=207;x2=238;
% y1=72;y2=97;
% hold on;
%
% plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1])
% title(['x1:',num2str(x1),', x2:',num2str(x2),', y1:',num2str(y1), ...
%     ', y2:',num2str(y2)])
%
x_saver2=x_saver2-x1+1;
y_saver2=y_saver2-y1+1;
 ROI1=img(y1:y2,x1:x2);
[cx,cy,cI]=improfile(ROI1,x_saver2,y_saver2);
bw=roipoly(ROI1,cx,cy);
%% set threshold
% single threshold is being used;
% double threshold method could be used in future for improvement
% threshold=.4;
% minROI=min(ROI1(:));
% maxROI=max(ROI1(:));
% threshold_true=minROI+(maxROI-minROI)*threshold;
%% whether to remove the single bright spot that could be noise?
singleRemove=1; % 1: remove single noise point; 0: not operation
searchRadius=1;
positiveNO=3;   % if the positive number is smaller than positiveNO, then remove it
%% whether to weight spots with corresponding grayscale value?
% scaleWeight=1;
% scaleWeight=0;
%% take a look at bright area
inputs.singleRemove=singleRemove;
inputs.searchRadius=searchRadius;
inputs.positiveNO=positiveNO;
% ROI_binary=imgBinary(ROI1,threshold_true,inputs);
%
% h_fig=figure;
% imagePlay(h_fig,BitDepth,ROI1,ROI_binary,threshold_true)
resultPath1=resultPath0;
result=[num2str(cellNO),'x1_',num2str(x1),'x2_',num2str(x2),'y1_',num2str(y1),'y2_',num2str(y2)];
resultPath=fullfile(resultPath1,result);
cd(resultPath1)
if exist(resultPath)==7
else
    mkdir(result)
end
%% get centroids
p=length(fileName);
centroids=zeros(p,2);
energy=zeros(p,1);
[m1,n1]=size(ROI1);
yyPositionTmp=(1:m1).';
yyPosition=yyPositionTmp*ones(1,n1);
xxPositionTmp=(1:n1).';
xxPosition=xxPositionTmp*ones(1,m1);
xxPosition=xxPosition.';
area=zeros(p,1);
for ii=1:p
    cd(filePath)
    img=differentTypeRead(fileName(ii).name,fileType);
    ROI1=img(y1:y2,x1:x2);
    % true threshold
    minROI=min(ROI1(:));
    maxROI=max(ROI1(:));
    threshold_true=minROI+(maxROI-minROI)*threshold;
    %
    ROI_binary=imgBinary(ROI1,threshold_true,inputs,othersH);
    ROI_binary=and(ROI_binary,bw);
    area(ii)=sum(ROI_binary(:));
    if scaleWeight~=1
        sumROI=sum(ROI_binary(:));
        yCentroidMat=yyPosition.*ROI_binary;
        xCentroidMat=xxPosition.*ROI_binary;
    else
        ROI_binary_weight=ROI_binary.*ROI1;
        sumROI=sum(ROI_binary_weight(:));
        yCentroidMat=yyPosition.*ROI_binary_weight;
        xCentroidMat=xxPosition.*ROI_binary_weight;
    end
    centroids(ii,1)=sum(yCentroidMat(:))/sumROI; % centroid at y axis
    centroids(ii,2)=sum(xCentroidMat(:))/sumROI; %centroid at x axis
    energy(ii)=sumROI;
    %% plot first 30 images
    if ii>0 && ii<=30
    h_fig=11;
    figure(h_fig);
    imagePlay(h_fig,BitDepth,ROI1,ROI_binary,threshold_true,centroids(ii,:))
    cd(resultPath)
    saveas(h_fig,num2str(ii,'%04d'),'tiff')
    end
   % close(h_fig)
end
cd(resultPath)
shift=centroids-ones(p,1)*centroids(1,:);
centroids_save=zeros(p,5);
centroids_save(:,1)=[1:p]';
centroids_save(:,2:3)=shift;
centroids2=shift;
centroids2(:,1)=IOS_time_gui_filter_shift(centroids2(:,1),filterMethod);
centroids2(:,2)=IOS_time_gui_filter_shift(centroids2(:,2),filterMethod);
centroids_save(:,4:5)=centroids2;
save([result,'_centroids.txt'],'centroids_save','-ASCII');
figure(4); subplot(2,1,1); plot(1:p,shift(:,2),'+',1:p,centroids2(:,2)); xlabel('horizontal shift')
subplot(2,1,2); plot(1:p,shift(:,1),'+',1:p,centroids2(:,1),'-'); xlabel('vertical shift');
saveas(gcf,'shift','tiff')
cd(filePath)
if nargout>=2
    varargout{1}=energy;
end
if nargout>=3
    varargout{2}=area;
end
function ROI_binary=imgBinary(ROI1,threshold_true,inputs,othersH)
if othersH.bigger==1
    ROI_binary=ROI1>=threshold_true;
else
    ROI_binary=ROI1<=threshold_true;
end
singleRemove=inputs.singleRemove;
if singleRemove==1 % remove single noise point
    positiveNO=inputs.positiveNO;
    searchRadius=inputs.searchRadius;
    [m,n]=size(ROI_binary);
    ROI_binary_tmp=false(m+2*searchRadius,n+2*searchRadius);
    ROI_binary_tmp(searchRadius+1:searchRadius+m,searchRadius+1:searchRadius+n)=ROI_binary;
    for ii=1:m
        for jj=1:n
            if ROI_binary(ii,jj)==1
            ROI_tmp=ROI_binary_tmp(ii:ii+2*searchRadius,jj:jj+2*searchRadius);
            ROI_binary(ii,jj)=sum(ROI_tmp(:))>=positiveNO;
            end
        end
    end
end
function img2=imageShow_centroid(h_axis,BitDepth,img)
axes(h_axis)
if BitDepth==8
    img2=uint8(img);
    imshow(img2)
% elseif BitDepth==16
%     img2=uint16(img);
%     imshow(img2)
else
    imgScale=1;
    img2=uint16(img*imgScale);
    imshow(img2);
    %imsho(im2uint16(mat2gray(imgScale)))
end
function imagePlay(h_fig,BitDepth,ROI1,ROI_binary,threshold_true,varargin)
figure(h_fig);
subplot(1,3,1);
ROI1_color_tmp=imageShow_centroid(gca,BitDepth,ROI1);
ROI1_color=zeros(size(ROI1,1),size(ROI1,2),3,'uint8');
ROI1_color(:,:,1)=ROI1_color_tmp;
ROI1_color(:,:,2)=ROI1_color_tmp;
ROI1_color(:,:,3)=ROI1_color_tmp;
ROI1_color_tmp(ROI_binary)=0;
ROI1_color(:,:,2)=ROI1_color_tmp;
ROI1_color(:,:,3)=ROI1_color_tmp;
subplot(1,3,2); imshow(ROI1_color);
hold on;
if nargin>5
    displacement=varargin{1};
    y=displacement(1);
    x=displacement(2);
plot(x,y,'g+','linewidth',2);
end
hold off;
subplot(1,3,3); hist(ROI1(:),20);
yRange=get(gca,'ylim');
hold on; plot(threshold_true*[1,1].',yRange.');
hold off;
function t1=IOS_time_gui_filter_shift(t1,filterMethod,varargin)
sizNo=filterMethod.sizNo;
sigmaNo=filterMethod.sigmaNo;
nameNo=filterMethod.name;
siz=sizNo;
sigma=sigmaNo;
width=siz;
[m1,n1]=size(t1);
if nargin==2
    if n1>m1
        t1=t1.';
    end
    if nameNo==0
        %imgAvg=imgAvg;
    elseif nameNo==1
        w=fspecial('gaussian',[siz,1],sigma);
        t1=imfilter(t1,w,'replicate');
    elseif nameNo==2
        t1=IOS_time_medianFilterAdapt(t1,siz);
    elseif nameNo==3
            t1=wiener2(t1,[siz,1]);
    elseif nameNo==4
        t1=edge(t1,'canny',[.05,.4],sigma);
    end
    if n1>m1
        t1=t1.';
    end
elseif nargin==3
    if nameNo==0
        %imgAvg=imgAvg;
    elseif nameNo==1
        w=fspecial('gaussian',[siz,1],sigma);
        t1=imfilter(t1,w,'replicate');
    elseif nameNo==2
        t1=medianFilterAdapt(t1,siz);
    elseif nameNo==3
            t1=wiener2(t1,[siz,1]);
    elseif nameNo==4
        t1=edge(t1,'canny',[.05,.4],sigma);
    end
end
function et_shiftEst_wt_Callback(hObject, eventdata, handles)
% hObject    handle to et_shiftEst_wt (see GCBO)
% et_shiftEst_wt as text
%        str2double(get(hObject,'String')) returns contents of et_shiftEst_wt as a double
function et_shiftEst_wt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_shiftEst_wt (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_cm.
function pb_cm_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cm (see GCBO)
file1=handles.file;
filePath1=handles.filePath;
cd(filePath1)
[file2,filePath2]=uigetfile('*.*','select query');
sizeYX=eval(get(handles.et_cm_size,'string'));
sizeY=sizeYX(1);
sizeX=sizeYX(2);
margin=sizeYX(3);
t1=cputime;
[yMap,xMap,ccMap]=sliceShift(file1,filePath1,file2,filePath2,sizeY,sizeX,margin,handles);
t2=cputime;
disp(['consuming time',num2str(t2-t1)])
function [yMap,xMap,ccMap]=sliceShift(file1,filePath1,file2,filePath2,sizeY,sizeX,margin,handles)
% [file1,filePath1]=uigetfile('*.*','select reference');
% cd(filePath1)
% [file2,filePath2]=uigetfile('*.*','select query');
img1=imread(fullfile(filePath1,file1));
img2=imread(fullfile(filePath2,file2));
% img1=img1(1:130,:);
% img2=img2(1:130,:);
% margin=3;
marginY=margin;
marginX=margin;
% sizeY=25;%% odds
% sizeX=25;%% odds
% sizeY=13;%% odds
% sizeX=7;%% odds
% sizeY=9;%% odds
% sizeX=7;%% odds
% sizeY=10;%% odds
% sizeX=10;%% odds
sizeY1=round((sizeY-1)/2);
sizeX1=round((sizeX-1)/2);
sizeY2=sizeY1+marginY;
sizeX2=sizeX1+marginX;
edgeY=marginY+sizeY;
edgeX=margin+sizeX;
[mm,nn]=size(img2);
%% make the images bigger; replicate padding
% img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
% img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
img1_big(1:edgeY,:)=ones(edgeY,1)*img1_big(1+edgeY,:);
img1_big(edgeY+mm+1:end,:)=ones(edgeY,1)*img1_big(mm,:);
img1_big(:,1:edgeX)=img1_big(:,edgeX+1)*ones(1,edgeX);
img1_big(:,nn+1+edgeX:end)=img1_big(:,nn)*ones(1,edgeX);
%% img2
img2_big=zeros(mm+2*edgeY,nn+2*edgeX);
img2_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img2;
img2_big(1:edgeY,:)=ones(edgeY,1)*img2_big(1+edgeY,:);
img2_big(edgeY+mm+1:end,:)=ones(edgeY,1)*img2_big(mm,:);
img2_big(:,1:edgeX)=img2_big(:,edgeX+1)*ones(1,edgeX);
img2_big(:,edgeX+nn+1:end)=img2_big(:,nn)*ones(1,edgeX);
%% initialize matrix
yMap=zeros(mm,nn);
xMap=zeros(mm,nn);
ccMap=zeros(mm,nn);
% absMap=zeros(mm,nn);
% phaseMap=zeros(mm,nn);
%% test peak
iiTmp=1;
ii=edgeY+iiTmp;
jjTmp=1;
jj=edgeX+jjTmp;
imgTrash=rand(size(img1_big));
ROI=imgTrash(ii-sizeY1:ii-sizeY1+sizeY-1,jj-sizeX1:jj-sizeX1+sizeX-1);
ROIplus=imgTrash(ii-sizeY2:ii+sizeY2,jj-sizeX2:jj+sizeX2);
[yPeak0,xPeak0]=shiftFinder(ROI,ROIplus);
%%
result1=['map_y',num2str(sizeY),'_x',num2str(sizeX),'_margin',num2str(margin)];
if get(handles.cb_cm_fast,'value')
    result1=['fast_',result1];
end
resultPath1=fullfile(filePath1,result1);
cd(filePath1)
if exist(resultPath1)==7 %#ok<*EXIST>
else
    mkdir(result1)
end
cd(resultPath1)
[fileTmptrash,resultPath1]=uiputfile('doNotChangeMyName.mat'); %#ok<*ASGLU>
h_wait=waitbar(0,'wait');
matlabpool
if get(handles.cb_cm_fast,'value')
    parfor iiTmp=1:mm
%          waitbar(iiTmp/mm,h_wait,[num2str(100*iiTmp/mm,'%04.1f'),'%completed']);
        ii=edgeY+iiTmp;
        for jjTmp=1:nn
            jj=edgeX+jjTmp;
            ROI=img1_big(ii-sizeY1:ii-sizeY1+sizeY-1,jj-sizeX1:jj-sizeX1+sizeX-1); %#ok<*PFBNS>
            ROI2=img2_big(ii-sizeY1:ii-sizeY1+sizeY-1,jj-sizeX1:jj-sizeX1+sizeX-1);
%             cc2=corrcoef(ROI(:),ROI2(:));
%             ccMap(iiTmp,jjTmp)=cc2(1,2);
            output= dftregistration(fft2(ROI),fft2(ROI2),margin);
            yMap(iiTmp,jjTmp)=-output(3);
            xMap(iiTmp,jjTmp)=-output(4);
    %         disp(['y',num2str(iiTmp,'%03d'),'jj',num2str(jjTmp,'%03d')])
        end
    end
else
    parfor iiTmp=1:mm
%          waitbar(iiTmp/mm,h_wait,[num2str(100*iiTmp/mm,'%04.1f'),'%completed']);
        ii=edgeY+iiTmp;
        for jjTmp=1:nn
            jj=edgeX+jjTmp;
            ROI=img1_big(ii-sizeY1:ii-sizeY1+sizeY-1,jj-sizeX1:jj-sizeX1+sizeX-1);
            ROIplus=img2_big(ii-sizeY2:ii+sizeY2,jj-sizeX2:jj+sizeX2);
            ROI2=img2_big(ii-sizeY1:ii-sizeY1+sizeY-1,jj-sizeX1:jj-sizeX1+sizeX-1);
%             cc2=corrcoef(ROI(:),ROI2(:));
%             ccMap(iiTmp,jjTmp)=cc2(1,2);
            if std(ROI(:))<0.001
                yPeak=yPeak0;
                xPeak=xPeak0;
            else
                [yPeak,xPeak]=shiftFinder(ROI,ROIplus);
            end
            yMap(iiTmp,jjTmp)=yPeak-yPeak0;
            xMap(iiTmp,jjTmp)=xPeak-xPeak0;
    %         disp(['y',num2str(iiTmp,'%03d'),'jj',num2str(jjTmp,'%03d')])
        end
    end
end
close(h_wait)
matlabpool close
cd(resultPath1)
save('yMap','yMap')
save('xMap','xMap')
% save('z ccMap','ccMap')
shift=xMap+1i*yMap;
amplitude=abs(shift); %#ok<*NASGU>
save('amplitude','amplitude')
phase=angle(shift)/pi*180;
save('phase','phase')
function [yPeak,xPeak]=shiftFinder(ROI,ROIplus)
 XCn2=normxcorr2(ROI, ROIplus);
 [ROIsize1,ROIsize2]=size(ROI);
 XCn=XCn2(ROIsize1:end-ROIsize1+1, ROIsize2:end-ROIsize2+1);
[max_cc, imax]=max(XCn(:));
[yPeak,xPeak]=ind2sub(size(XCn),imax(1));
function [output Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007
% Portions of this code were taken from code written by Ann M. Kowalczyk
% and James R. Fienup.
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
% (1990).
% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).
% Inputs
% buf1ft    Fourier transform of reference image,
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register,
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.
% Default usfac to 1
if exist('usfac')~=1, usfac=1; end
% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft)));
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    output=[error,diffphase];
% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc);
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    md2 = fix(m/2);
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
% Partial-pixel shift
else
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
    % Compute crosscorrelation and locate the peak
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac;
        col_shift = round(col_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid
        [max1,loc1] = max(CC);
        [max2,loc2] = max(max1);
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
    % If upsampling = 2, no additional pixel shift refinement
    else
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end
% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
return
function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1)
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06
% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the
%     [roff+1 coff+1] element.
% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]
[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return
function et_cm_size_Callback(hObject, eventdata, handles)
function et_cm_size_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_eval_Callback(hObject, eventdata, handles)
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles)
function et_eval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in cb_eval.
function cb_eval_Callback(hObject, eventdata, handles)
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
function pb_save_Callback(hObject, eventdata, handles)
cd(handles.filePath)
handles.newIntensity=inline(get(handles.et_inline,'string'));
fileType=handles.fileType;
file=handles.file;
x=differentTypeRead(file,fileType);
d3=size(x,3);
fileName=handles.fileName;
p=length(fileName);
filePath=handles.filePath;
result1='eval';
resultPath1=fullfile(filePath,result1);
cd(filePath)
if exist(resultPath1)==7
else
    mkdir(result1)
end
cd(resultPath1)
% [fileTmptrash,resultPath1]=uiputfile(file);
[fileTmptrash,resultPath1,filterIndex]=uiputfile({'*.tif';'*.mat';'*.png'},'save as',file);

h_wait=waitbar(0,'wait');
if filterIndex==2
    for ii=1:p
        waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
    %     cd(filePath)
        file1=fileName(ii).name;
        x=differentTypeRead(fullfile(filePath,file1),fileType);
        x=handles.newIntensity(x);
        if get(handles.cb_eval,'value')
             eval(get(handles.et_eval,'string')) ;
        end
    %     cd(resultPath1)
        if strcmp(fileType,'mat')==0
              file1=[file1(1:end-length(fileType)),'mat'];
        end
        save(fullfile(resultPath1,file1),'x')
        
%             if get(handles.cb_uint16,'value')
%                 imwrite(uint16(x),fullfile(resultPath1,file1))
% 
%             else
%                 imwrite(uint8(x),fullfile(resultPath1,file1))
% 
%             end        

    end       
elseif filterIndex==1
    for ii=1:p
        waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
    %     cd(filePath)
        file1=fileName(ii).name;
        x=differentTypeRead(fullfile(filePath,file1),fileType);
        x=handles.newIntensity(x);
        if get(handles.cb_eval,'value')
             eval(get(handles.et_eval,'string')) ;
        end
    %     cd(resultPath1)
        if strcmp(fileType,'mat')
              file1=[file1(1:end-length(fileType)),'tif'];
        end
            if get(handles.cb_uint16,'value')
                if strcmp(file1(end-2:end),'tif')
                     imwrite(uint16(x),fullfile(resultPath1,file1),'compression','none')
                else
                    imwrite(uint16(x),fullfile(resultPath1,file1))
                end
                

            else
                if strcmp(file1(end-2:end),'tif')
                     imwrite(uint8(x),fullfile(resultPath1,file1),'compression','none')
                else
                    imwrite(uint8(x),fullfile(resultPath1,file1))
                end

            end        

    end   
elseif filterIndex==3
    for ii=1:p
        waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
    %     cd(filePath)
        file1=fileName(ii).name;
        x=differentTypeRead(fullfile(filePath,file1),fileType);
        x=handles.newIntensity(x);
        if get(handles.cb_eval,'value')
             eval(get(handles.et_eval,'string')) ;
        end
    %     cd(resultPath1)
        if strcmp(fileType,'mat')||strcmp(fileType,'tif')
              file1=[file1(1:end-length(fileType)),'png'];
        end
            if get(handles.cb_uint16,'value')
                 imwrite(uint16(x),fullfile(resultPath1,file1))

                

            else
                imwrite(uint8(x),fullfile(resultPath1,file1))
               
            end        

    end       
end
  close(h_wait)
% --- Executes on button press in pb_g7.
function pb_g7_Callback(hObject, eventdata, handles)
% hObject    handle to pb_g7 (see GCBO)
gN=7;
handles=getPointG1(handles,gN);
guidata(hObject,handles);
function pb_g7_et_Callback(hObject, eventdata, handles)
function pb_g7_et_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_phase_sect_Callback(hObject, eventdata, handles)
function et_phase_sect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_phase_run.
function pb_phase_run_Callback(hObject, eventdata, handles)
% hObject    handle to pb_phase_run (see GCBO)
uiwait(msgbox('make sure 1st is xMap;2nd YMap','Title','modal'));
y0x0=eval(get(handles.et_phase_y0x0,'string'));
y0=y0x0(1);y0=round(y0);
x0=y0x0(2);x0=round(x0);
sectionTmp=eval(get(handles.et_phase_sect,'string'));
section=sectionTmp(1);
initial=sectionTmp(2);
lengthAmp=sectionTmp(3);
filePath=handles.filePath;
fileType=handles.fileType;
file=handles.file;
fileName=handles.fileName;
cd(filePath)
%% creating secitons
img=differentTypeRead(file,fileType);
[m1,n1,o1]=size(img);
rMax=floor(min([y0,m1-y0+1,x0,n1-x0+1]));
theta=linspace(0,2*pi,section+1)+initial/180*pi;
% thetaY=y0+cos(theta)*rMax;
% thetaX=x0+sin(theta)*rMax;
sectionXY=cell(section,1);
step=0.5/180*pi;
% edgeStep=step;
edgeStep=0;
pointNO=round(2*pi/step/section);
   handles.newIntensity=inline(get(handles.et_inline,'string'));
    newImg=handles.newIntensity(img);
    newImg=uint8(newImg);
    figure(90);close(90);figure(90);imshow(newImg);hold on;
for ii=1:section
    thetaTmp=linspace(theta(ii),theta(ii+1)-edgeStep,pointNO);
    thetaY=y0+sin(thetaTmp)*rMax;
    thetaX=x0+cos(thetaTmp)*rMax;
    sectionXY{ii}{1}=[x0,thetaX,x0].';
    sectionXY{ii}{2}=[y0,thetaY,y0].';
    plot(sectionXY{ii}{1},sectionXY{ii}{2},'r','lineWidth',2)
end
%% imread xMap, yMap
xMap=differentTypeRead(fileName(1).name,fileType);
yMap=differentTypeRead(fileName(2).name,fileType);
%%
xMap2=xMap;
yMap2=yMap;
xMap(abs(xMap2)>5)=0;
yMap(abs(yMap2)>5)=0;
xMap(abs(yMap2)>5)=0;
yMap(abs(xMap2)>5)=0;
handles.phaseH.xMap=xMap;
handles.phaseH.yMap=yMap;
%% normalization vectors;
amplitude=abs(xMap+1i*yMap);
ampLogic=logical(amplitude);
xMap(ampLogic)=xMap(ampLogic)./(amplitude(ampLogic));
yMap(ampLogic)=yMap(ampLogic)./(amplitude(ampLogic));
%%
pointPst=handles.pointPst;
cx=pointPst{7}{1}(:,1);
cy=pointPst{7}{1}(:,2);
bw=roipoly(xMap(:,:,1),cx,cy);
plot(cx,cy,'g','lineWidth',2)
% xMap=xMap.*bw;
% yMap=yMap.*bw;
xyMap=xMap+1i*yMap;
phaseMap=zeros(section,1);
phaseMapX=zeros(section,1);
phaseMapY=zeros(section,1);
%%
yyPositionTmp=(1:m1).';
yyPosition=yyPositionTmp*ones(1,n1);
xxPositionTmp=(1:n1).';
xxPosition=xxPositionTmp*ones(1,m1);
xxPosition=xxPosition.';
for ii=1:section
    cx2=sectionXY{ii}{1};
    cy2=sectionXY{ii}{2};
    bw2=roipoly(xMap(:,:,1),cx2,cy2);
    bwTmp=and(bw,bw2);
    sumBW=sum(bwTmp(:));
    centroidY=yyPosition.*bwTmp;
    centroidX=xxPosition.*bwTmp;
    phaseMapY(ii)=sum(centroidY(:))/sumBW;
    phaseMapX(ii)=sum(centroidX(:))/sumBW;
    plot(phaseMapX(ii),phaseMapY(ii),'r+','lineWidth',2)
    bw3=and(bwTmp,ampLogic);
    xy=xyMap(bw3);
    phaseMap(ii)=sum(xy);
    %%
    if sum(bw3(:))<=3
      phaseMap(ii)=0;
    end
    %%
    if abs(phaseMap(ii))~=0
        phaseMap(ii)=phaseMap(ii)/abs(phaseMap(ii));
         quiver(phaseMapX(ii),phaseMapY(ii),real(phaseMap(ii)),imag(phaseMap(ii)),lengthAmp)
    end
%     for jj=1:length(xy)
%         phaseMap(ii)=phaseMap(ii)+xy(jj);
%
%         phaseMap(ii)=phaseMap(ii)/abs(phaseMap(ii));
%     end
end
%  quiver(phaseMapX,phaseMapY,real(phaseMap),imag(phaseMap),50)
handles.phaseH.phaseMap=phaseMap;
handles.phaseH.phaseMapX=phaseMapX;
handles.phaseH.phaseMapY=phaseMapY;
handles.phaseH.sectionXY=sectionXY;
handles.phaseH.cx=cx;
handles.phaseH.cy=cy;
handles.phaseH.x0=x0;
handles.phaseH.y0=y0;
%%
[m1,n1,o1]=size(img);
r1r2=eval(get(handles.et_phase_r1r2,'string'));
r1=r1r2(1);
r2=r1r2(2);
r3=r1r2(3);
% rMax=floor(min([y0,m1-y0+1,x0,n1-x0+1]));
theta=linspace(0,2*pi,section+1)+initial/180*pi;
% thetaY=y0+cos(theta)*rMax;
% thetaX=x0+sin(theta)*rMax;
sectionXY=cell(section,1);
    newImg=handles.newIntensity(img);
    newImg=uint8(newImg);
    figure(91);close(91);figure(91);imshow(newImg);hold on;
    thetaCenter=theta(1:end-1)+theta(2:end);
    thetaCenter=thetaCenter/2;
    phaseMapX=zeros(section,1);
    phaseMapY=zeros(section,2);
for ii=1:section
    y3=y0+sin(thetaCenter(ii))*r3;phaseMapY(ii)=y3;
    x3=x0+cos(thetaCenter(ii))*r3;phaseMapX(ii)=x3;
    thetaTmp=linspace(theta(ii),theta(ii+1)-edgeStep,pointNO);
    thetaY=y0+sin(thetaTmp)*r1;
    thetaX=x0+cos(thetaTmp)*r1;
    thetaY2=y0+sin(thetaTmp)*r2;
    thetaX2=x0+cos(thetaTmp)*r2;
    sectionXY{ii}{1}=[thetaX,thetaX2(end:-1:1),thetaX(1)].';
    sectionXY{ii}{2}=[thetaY,thetaY2(end:-1:1),thetaY(1)].';
    plot(sectionXY{ii}{1},sectionXY{ii}{2},'r','lineWidth',2)
    plot(x3,y3,'r+','lineWidth',1)
    if abs(phaseMap(ii))~=0
%         phaseMap(ii)=phaseMap(ii)/abs(phaseMap(ii));
         quiverH=quiver(phaseMapX(ii),phaseMapY(ii),1*real(phaseMap(ii)),1*imag(phaseMap(ii)),lengthAmp);
         set(quiverH,'lineWidth',1,'color',[0,1,0])
         adjust_quiver_arrowhead_size(quiverH,5);
%          set(quiverH,'autoScale','off')
    end
end
% figure(92);close(92);figure(92);
 img92=zeros(m1,n1);%imshow(uint8(img92)); hold on;
%%
for ii=1:section
    bw=roipoly(img,sectionXY{ii}{1},sectionXY{ii}{2});
%     plot(sectionXY{ii}{1},sectionXY{ii}{2},'r','lineWidth',2)
  %  plot(x3,y3,'r+','lineWidth',1)
    if abs(phaseMap(ii))~=0
%         phaseMap(ii)=phaseMap(ii)/abs(phaseMap(ii));
%          quiverH=quiver(phaseMapX(ii),phaseMapY(ii),1*real(phaseMap(ii)),1*imag(phaseMap(ii)),lengthAmp);
%          set(quiverH,'lineWidth',1,'color',[0,1,0])
%          adjust_quiver_arrowhead_size(quiverH,5);
         ampTemp=amplitude(and(ampLogic,bw));
         img92(bw)=mean(ampTemp);
%          set(quiverH,'autoScale','off')
    end
end
figure(93);imshow(im2uint8(mat2gray(img92,[-3,3])));set(gcf,'colormap',jet(126))
handles.phaseH.sectionXY2=sectionXY;
handles.phaseH.phaseMapX2=phaseMapX;
handles.phaseH.phaseMapY2=phaseMapY;
handles.phaseH.ampMap=img92;
guidata(hObject,handles)
function adjust_quiver_arrowhead_size(quivergroup_handle, scaling_factor)
% Make quiver arrowheads bigger or smaller.
%
% adjust_quiver_arrowhead_size(quivergroup_handle, scaling_factor)
%
% Example:
%   h = quiver(1:100, 1:100, randn(100, 100), randn(100, 100));
%   adjust_quiver_arrowhead_size(h, 1.5);   % Makes all arrowheads 50% bigger.
%
% Inputs:
%   quivergroup_handle      Handle returned by "quiver" command.
%   scaling_factor          Factor by which to shrink/grow arrowheads.
%
% Output: none
% Kevin J. Delaney
% December 21, 2011
% BMT Scientific Marine Services (www.scimar.com)
if ~exist('quivergroup_handle', 'var')
    help(mfilename);
    return
end
if isempty(quivergroup_handle) || any(~ishandle(quivergroup_handle))
    errordlg('Input "quivergroup_handle" is empty or contains invalid handles.', ...
             mfilename);
    return
end
if length(quivergroup_handle) > 1
    errordlg('Expected "quivergroup_handle" to be a single handle.', mfilename);
    return
end
if ~strcmpi(get(quivergroup_handle, 'Type'), 'hggroup')
    errrodlg('Input "quivergroup_handle" is not of type "hggroup".', mfilename);
    return
end
if ~exist('scaling_factor', 'var') || ...
   isempty(scaling_factor) || ...
   ~isnumeric(scaling_factor)
    errordlg('Input "scaling_factor" is missing, empty or non-numeric.', ...
             mfilename);
    return
end
if length(scaling_factor) > 1
    errordlg('Expected "scaling_factor" to be a scalar.', mfilename);
    return
end
if scaling_factor <= 0
    errordlg('"Scaling_factor" should be > 0.', mfilename);
    return
end
line_handles = get(quivergroup_handle, 'Children');
if isempty(line_handles) || (length(line_handles) < 3) || ...
   ~ishandle(line_handles(2)) || ~strcmpi(get(line_handles(2), 'Type'), 'line')
    errordlg('Unable to adjust arrowheads.', mfilename);
    return
end
arrowhead_line = line_handles(2);
XData = get(arrowhead_line, 'XData');
YData = get(arrowhead_line, 'YData');
if isempty(XData) || isempty(YData)
    return
end
%   Break up XData, YData into triplets separated by NaNs.
first_nan_index = find(~isnan(XData), 1, 'first');
last_nan_index  = find(~isnan(XData), 1, 'last');
for index = first_nan_index : 4 : last_nan_index
    these_indices = index + (0:2);
    if these_indices(end) > length(XData)
        break
    end
    x_triplet = XData(these_indices);
    y_triplet = YData(these_indices);
    if any(isnan(x_triplet)) || any(isnan(y_triplet))
        continue
    end
    %   First pair.
    delta_x = diff(x_triplet(1:2));
    delta_y = diff(y_triplet(1:2));
    x_triplet(1) = x_triplet(2) - (delta_x * scaling_factor);
    y_triplet(1) = y_triplet(2) - (delta_y * scaling_factor);
    %   Second pair.
    delta_x = diff(x_triplet(2:3));
    delta_y = diff(y_triplet(2:3));
    x_triplet(3) = x_triplet(2) + (delta_x * scaling_factor);
    y_triplet(3) = y_triplet(2) + (delta_y * scaling_factor);
    XData(these_indices) = x_triplet;
    YData(these_indices) = y_triplet;
end
set(arrowhead_line, 'XData', XData, 'YData', YData);
function et_phase_y0x0_Callback(hObject, eventdata, handles)
function et_phase_y0x0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pb_phase_run_CreateFcn(hObject, eventdata, handles)
function pb_phase_show_Callback(hObject, eventdata, handles)
% hObject    handle to pb_phase_show (see GCBO)
phaseH=handles.phaseH;
sectionXY=phaseH.sectionXY2;
phaseMapX=phaseH.phaseMapX2;
phaseMapY=phaseH.phaseMapY2;
phaseMap=phaseH.phaseMap;
x0=phaseH.x0;
y0=phaseH.y0;
section=length(sectionXY);
filePath=handles.filePath;
fileType=handles.fileType;
file=handles.file;
fileName=handles.fileName;
cd(filePath)
sectionTmp=eval(get(handles.et_phase_sect,'string'));
% section=sectionTmp(1);
% initial=sectionTmp(2);
lengthAmp=sectionTmp(3);
%% creating secitons
img=differentTypeRead(file,fileType);
[m1,n1,o1]=size(img);
   handles.newIntensity=inline(get(handles.et_inline,'string'));
    newImg=handles.newIntensity(img);
    newImg=uint8(newImg);
    figure(91);close(91);figure(91);imshow(newImg);hold on;
for ii=1:section
    plot(sectionXY{ii}{1},sectionXY{ii}{2},'r','lineWidth',2)
    plot(phaseMapX(ii),phaseMapY(ii),'r+','lineWidth',1)
    if abs(phaseMap(ii))~=0
         quiverH=quiver(phaseMapX(ii),phaseMapY(ii),1*real(phaseMap(ii)),1*imag(phaseMap(ii)),lengthAmp);
         set(quiverH,'lineWidth',1,'color',[0,1,0])
         adjust_quiver_arrowhead_size(quiverH,5);
    end
end
result1=['phase'];
resultPath1=fullfile(filePath,result1);
cd(filePath)
if exist(resultPath1)==7
else
    mkdir(result1)
end
cd(resultPath1)
[file,resultPath1]=uiputfile('phaseHandle');
cd(resultPath1)
save(file,'phaseH');
saveas(gcf,file,'tiff')
function et_phase_r1r2_Callback(hObject, eventdata, handles)
function et_phase_r1r2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton119_Callback(hObject, eventdata, handles)
uiwait(msgbox('make sure 1st is xMap;2nd YMap','Title','modal'));
fileName=handles.fileName;
filePath=handles.filePath;
file=handles.file;
p=length(fileName);
fileType=handles.fileType;
fileTypeLength=length(fileType);
%% binning size
xBinning=eval(get(handles.et_xBinning,'string'));
yBinning=eval(get(handles.et_yBinning,'string'));
%% frame rate unit HZ
frameRate=1000;
%% ROI
xStart=eval(get(handles.et_xBinningStart,'string'));
xEnd=eval(get(handles.et_xBinningEnd,'string'));
yStart=eval(get(handles.et_yBinningStart,'string'));
yEnd=eval(get(handles.et_yBinningEnd,'string'));
%% verify whether ROI fits
if mod(yEnd-yStart+1,yBinning)~=0
    yEnd=floor((yEnd-yStart+1)/yBinning)*yBinning-1+yStart;
end
if mod(xEnd-xStart+1,xBinning)~=0
    xEnd=floor((xEnd-xStart+1)/xBinning)*xBinning-1+xStart;
end
new_xsize=(xEnd-xStart+1)/xBinning; % image size after pixel binning
new_ysize=(yEnd-yStart+1)/yBinning; % image size after pixel binning
%% file name for result
name=['phaseRect_x',num2str(xBinning),'y',num2str(yBinning)];
quiverSize=eval(get(handles.et_quiverSize,'string'));
name=[name,'_',get(handles.et_xy_noise,'string'),'_',num2str(quiverSize(1)),'_',num2str(quiverSize(2))];
result0='phase';
IOSresultPath=get(handles.et_IOSresultPath,'string');
resultPath0=fullfile(IOSresultPath,result0);
cd(IOSresultPath)
if exist(resultPath0)==7
else
    mkdir(result0)
end
cd(resultPath0)
%% subdirectory
result=name;
resultPath=fullfile(resultPath0,result);
cd(resultPath0)
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
[fileTmp,resultPath]=uiputfile('phaseH.mat');
%%
cd(filePath)
img=differentTypeRead(handles.file,fileType);
[m1,n1,o1]=size(img);
handles.newIntensity=inline(get(handles.et_inline,'string'));

newImg=handles.newIntensity(img);
newImg=uint8(newImg);
figure(90);close(90);figure(90);imshow(newImg);hold on;
%% imread xMap, yMap
xMap=differentTypeRead(fileName(1).name,fileType);
yMap=differentTypeRead(fileName(2).name,fileType);
phaseH.xMap=xMap;
phaseH.yMap=yMap;
% xMap=xMap2(yStart:yEnd,xStart:xEnd);
% yMap=yMap2(yStart:yEnd,xStart:xEnd);
% handles.phaseH_slice.xMap=xMap;
% handles.phaseH_slice.yMap=yMap;
%% normalization vectors;
amplitude=abs(xMap+1i*yMap);
ampLogic=logical(amplitude);
xMap(ampLogic)=xMap(ampLogic)./(amplitude(ampLogic));
yMap(ampLogic)=yMap(ampLogic)./(amplitude(ampLogic));
xyMap=xMap+1i*yMap;
%%
xSect=1:xBinning:(1+new_xsize)*xBinning;
xSect=xSect+xStart-1;
xSect(end)=xSect(end)-1;
ySect=1:yBinning:(1+new_ysize)*yBinning;
ySect=ySect+yStart-1;
ySect(end)=ySect(end)-1;
cx=cell(new_ysize,new_xsize);
cy=cell(new_ysize,new_xsize);
phaseMapX=zeros(new_ysize,new_xsize);
phaseMapY=zeros(new_ysize,new_xsize);
phaseMap=zeros(new_ysize,new_xsize);
edgeStep=0;
for ii=1:new_ysize
    for jj=1:new_xsize
        cx{ii,jj}=[xSect(jj),xSect(jj),xSect(jj+1)-edgeStep,xSect(jj+1)-edgeStep,xSect(jj)].';
        cy{ii,jj}=[ySect(ii),ySect(ii+1)-edgeStep,ySect(ii+1)-edgeStep,ySect(ii),ySect(ii)].';
        phaseMapY(ii,jj)=0.5*(ySect(ii)+ySect(ii+1));
         phaseMapX(ii,jj)=0.5*(xSect(jj)+xSect(jj+1));
        plot(cx{ii,jj},cy{ii,jj},'r-','linewidth',1)
%         plot(phaseMapX(ii,jj),phaseMapY(ii,jj),'r+','linewidth',1)
    end
end
%%
yyPositionTmp=(1:m1).';
yyPosition=yyPositionTmp*ones(1,n1);
xxPositionTmp=(1:n1).';
xxPosition=xxPositionTmp*ones(1,m1);
xxPosition=xxPosition.';
mapAll=zeros(m1,n1);
for ii=1:new_ysize
    for jj=1:new_xsize
        cx2=cx{ii,jj};
        cy2=cy{ii,jj};
        bw2=roipoly(xMap(:,:,1),cx2,cy2);
%         bwTmp=and(bw,bw2);
        bwTmp=bw2;
        bwXY=and(bwTmp,ampLogic);
        xy=xyMap(bwXY);
        if sum(bwXY(:))>eval(get(handles.et_xy_noise,'string'))
            phaseMap(ii,jj)=sum(xy);
        end
    if abs(phaseMap(ii,jj))~=0
        phaseMap(ii,jj)=phaseMap(ii,jj)/abs(phaseMap(ii,jj));
        quiverH=quiver(phaseMapX(ii,jj),phaseMapY(ii,jj),real(phaseMap(ii,jj)),imag(phaseMap(ii,jj)),quiverSize(1));
        set(quiverH,'lineWidth',1,'color',[0,1,0])
        adjust_quiver_arrowhead_size(quiverH,quiverSize(2));
        mapAll(bw2)=mean(amplitude(bwXY));
    end
    end
end
mapAllTmp=mapAll(mapAll>0);
resizeFactor=1;
if resizeFactor==2
% upperLim=2*mean(mapAllTmp);
    upperLim=1;
    % figure(91);close(91);figure(91);imshow(im2uint8(mat2gray(mapAll,[-upperLim,upperLim])));hold on;
    mapAll=mapAll/2;
     if get(handels.cb_map_raw,'value')
         amplitude=amplitude/2;
        figure(91);close(91);figure(91);imshow(amplitude,[0,upperLim]);hold on;
    else
        figure(91);close(91);figure(91);imshow(mapAll,[0,upperLim]);hold on;
    end   
%     figure(91);close(91);figure(91);imshow(mapAll,[0,upperLim]);hold on;
    map0=jet(256*2);
    map=map0(256:end,:);
    set(gcf,'colormap',map)
    hh=colorbar;
    set(hh,'yTick',linspace(0,upperLim,4))
    % set(hh,'YTick',linspace(0,2,4))
    set(hh,'YTickLabel',{'0','0.1','0.2','0.3'})
    set(hh,'FontSize',15)
else
    upperLim=2;
    if get(handles.cb_map_raw,'value')
        figure(91);close(91);figure(91);imshow(amplitude,[0,upperLim]);hold on;
    else
        figure(91);close(91);figure(91);imshow(mapAll,[0,upperLim]);hold on;
    end
    % figure(91);close(91);figure(91);imshow(im2uint8(mat2gray(mapAll,[-upperLim,upperLim])));hold on;
%     mapAll=mapAll/2;
    
    map0=jet(256*2);
    map=map0(256:end,:);
    set(gcf,'colormap',map)
    hh=colorbar;
    set(hh,'yTick',linspace(0,upperLim,4))
    % set(hh,'YTick',linspace(0,2,4))
    set(hh,'YTickLabel',{'0','0.2','0.4','0.6'})
    set(hh,'FontSize',15)    
end
for ii=1:new_ysize
    for jj=1:new_xsize
        cx2=cx{ii,jj};
        cy2=cy{ii,jj};
        bw2=roipoly(xMap(:,:,1),cx2,cy2);
%         bwTmp=and(bw,bw2);
        bwTmp=bw2;
        bwXY=and(bwTmp,ampLogic);
        xy=xyMap(bwXY);
        if sum(bwXY(:))>eval(get(handles.et_xy_noise,'string'))
            phaseMap(ii,jj)=sum(xy);
        end
    if abs(phaseMap(ii,jj))~=0
        phaseMap(ii,jj)=phaseMap(ii,jj)/abs(phaseMap(ii,jj));
%         phaseMap(ii,jj)=phaseMap(ii,jj)/abs(phaseMap(ii,jj));
            theta=angle(phaseMap(ii,jj))+pi;
            phaseMapX(ii,jj)=phaseMapX(ii,jj)+cos(theta)*0.4*xBinning;
            phaseMapY(ii,jj)=phaseMapY(ii,jj)+sin(theta)*0.4*xBinning;           
        quiverH=quiver(phaseMapX(ii,jj),phaseMapY(ii,jj),real(phaseMap(ii,jj)),imag(phaseMap(ii,jj)),quiverSize(1));
        set(quiverH,'lineWidth',2.5,'color',[0,0,1])
        adjust_quiver_arrowhead_size(quiverH,quiverSize(2));
    end
    plot(cx{ii,jj},cy{ii,jj},'r-','linewidth',1)
    end
end
cd(resultPath)
phaseH.phaseMapX=phaseMapX;
phaseH.phaseMapY=phaseMapY;
phaseH.cx=cx;
phaseH.cy=cy;
phaseH.phaseMap=phaseMap;
phaseH.ampMap=mapAll;
save(fileTmp,'phaseH')
a=(get(handles.et_xy_noise,'string'));

% saveas(91,['amp_dir',a,num2str(quiverSize(1)),'_',num2str(quiverSize(2))],'tiff')
% saveas(90,['direction',a,num2str(quiverSize(1)),'_',num2str(quiverSize(2))],'tiff')
    rePath_export_figs;
    figure(91);
%     cd(resultPath)
    colorbar off;
    
       set(gca,'color','none')
   export_fig(91,[['amp_dir',a,num2str(quiverSize(1)),'_',num2str(quiverSize(2))],'.png'],'-transparent')
   
% saveas(90,'direction','tiff')
% --- Executes on button press in cb_sc.
function cb_sc_Callback(hObject, eventdata, handles)
% --- Executes on button press in cb_cm_fast.
function cb_cm_fast_Callback(hObject, eventdata, handles)
function pushbutton120_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton120 (see GCBO)
sizeYX=eval(get(handles.et_cm_size,'string'));
sizeY=sizeYX(1);
sizeX=sizeYX(2);
margin=sizeYX(3);
%%
fileName=handles.fileName;
fileType=handles.fileType;
p=length(fileName);
file=handles.file;
filePath=handles.filePath;
cd(filePath)
img=differentTypeRead(file,fileType);
[mm,nn]=size(img);
img1=zeros(mm,nn);
%%
marginY=margin;
marginX=margin;
sizeY1=round((sizeY-1)/2); 
sizeX1=round((sizeX-1)/2);
sizeY2=sizeY1+marginY;
sizeX2=sizeX1+marginX;
edgeY=marginY+sizeY;
edgeX=margin+sizeX;
%% reference
p2=eval(get(handles.et_pre2,'string'));
for ii=1:p2
    img=differentTypeRead(fileName(ii).name,fileType);
    img1=img1+img;
end
img1=img1/p2;
img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
img1_big(1:edgeY,:)=ones(edgeY,1)*img1_big(1+edgeY,:);
img1_big(edgeY+mm+1:end,:)=ones(edgeY,1)*img1_big(mm,:);
img1_big(:,1:edgeX)=img1_big(:,edgeX+1)*ones(1,edgeX);
img1_big(:,nn+1+edgeX:end)=img1_big(:,nn)*ones(1,edgeX);
figure;imshow(im2uint8(mat2gray(img1_big)))
%%
    %%
    result1=['map_y',num2str(sizeY),'_x',num2str(sizeX),'_margin',num2str(margin)];
    if get(handles.cb_cm_fast,'value')
        result1=['fast_',result1];
    end
    resultPath1=fullfile(filePath,result1);
    cd(filePath)
    if exist(resultPath1)==7
    else
        mkdir(result1)
    end
    %% XMAP FOLDER
    cd(resultPath1)
    resultX='xMap';
    resultPathX=fullfile(resultPath1,resultX);
    if exist(resultPathX)==7
    else
        mkdir(resultX)
    end
    %% YMAP FOLDER
    cd(resultPath1)
    resultY='yMap';
    resultPathY=fullfile(resultPath1,resultY);
    if exist(resultPathY)==7
    else
        mkdir(resultY)
    end
    %% CC MAP FOLDER
%     cd(resultPath1)
%     resultCC='CCMap';
%     resultPathCC=fullfile(resultPath1,resultCC);
%     if exist(resultPathCC)==7
%     else
%         mkdir(resultCC)
%     end
    %% AMP MAP FOLDER
    cd(resultPath1)
    resultAMP='AMPMap';
    resultPathAMP=fullfile(resultPath1,resultAMP);
    if exist(resultPathAMP)==7
    else
        mkdir(resultAMP)
    end
%     [fileTmptrash,resultPath1]=uiputfile('results.mat');
%%
 h_wait=waitbar(0,'wait');
 t1=cputime;
 P=p;
     %% test peak
iiTmp=1;
ii=edgeY+iiTmp;
jjTmp=1;
jj=edgeX+jjTmp;
imgTrash=rand(size(img1_big));
ROI=imgTrash(ii-sizeY1:ii+sizeY1,jj-sizeX1:jj+sizeX1);
ROIplus=imgTrash(ii-sizeY2:ii+sizeY2,jj-sizeX2:jj+sizeX2);
[yPeak0,xPeak0]=shiftFinder(ROI,ROIplus);
matlabpool
for ss=1:p
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%',datestr(now)])
    if mod(ss,10)==0
        display('let CPU rest for 3 munites')
        %pause(60*3)
    end
   waitbar(ss/P,h_wait,[num2str(100*ss/P,'%04.1f'),'%completed']);
    %% make the images bigger; replicate padding
    % img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
    % img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
    %% img2
    cd(filePath)
    img2=differentTypeRead(fileName(ss).name,fileType);
    img2_big=zeros(mm+2*edgeY,nn+2*edgeX);
    img2_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img2;
    img2_big(1:edgeY,:)=ones(edgeY,1)*img2_big(1+edgeY,:);
    img2_big(edgeY+mm+1:end,:)=ones(edgeY,1)*img2_big(mm,:);
    img2_big(:,1:edgeX)=img2_big(:,edgeX+1)*ones(1,edgeX);
    img2_big(:,edgeX+nn+1:end)=img2_big(:,nn)*ones(1,edgeX);
    %% initialize matrix
    yMap=zeros(mm,nn);
    xMap=zeros(mm,nn);
%     ccMap=zeros(mm,nn);
    % absMap=zeros(mm,nn);
    % phaseMap=zeros(mm,nn);
    if get(handles.cb_cm_fast,'value')
        for iiTmp=1:mm
            ii=edgeY+iiTmp;
            img1_slice=img1_big(ii-sizeY1:ii+sizeY1,:);
            img2_slice=img2_big(ii-sizeY1:ii+sizeY1,:);
            parfor jjTmp=1:nn
                jj=edgeX+jjTmp;
                ROI=img1_slice(:,jj-sizeX1:jj+sizeX1);
                ROI2=img2_slice(:,jj-sizeX1:jj+sizeX1);
%                 cc2=corrcoef(ROI(:),ROI2(:));
%                 ccMap(iiTmp,jjTmp)=cc2(1,2);
                output= dftregistration(fft2(ROI),fft2(ROI2),margin);
                yMap(iiTmp,jjTmp)=-output(3);
                xMap(iiTmp,jjTmp)=-output(4);
        %         disp(['y',num2str(iiTmp,'%03d'),'jj',num2str(jjTmp,'%03d')])
            end
        end
    else
        for iiTmp=1:mm
            ii=edgeY+iiTmp;
            img1_slice=img1_big(ii-sizeY1:ii+sizeY1,:);
            img2_slice=img2_big(ii-sizeY2:ii+sizeY2,:);
            img2_slice2=img2_big(ii-sizeY1:ii+sizeY1,:);
            parfor jjTmp=1:nn
                jj=edgeX+jjTmp;
                ROI=img1_slice(:,jj-sizeX1:jj+sizeX1);
                ROIplus=img2_slice(:,jj-sizeX2:jj+sizeX2);
                ROI2=img2_slice2(:,jj-sizeX1:jj+sizeX1);
%                 cc2=corrcoef(ROI(:),ROI2(:));
%                 ccMap(iiTmp,jjTmp)=cc2(1,2);
                if std(ROI(:))<0.001
                    yPeak=yPeak0;
                    xPeak=xPeak0;
                else
                    [yPeak,xPeak]=shiftFinder(ROI,ROIplus);
                end
                yMap(iiTmp,jjTmp)=yPeak-yPeak0;
                xMap(iiTmp,jjTmp)=xPeak-xPeak0;
        %         disp(['y',num2str(iiTmp,'%03d'),'jj',num2str(jjTmp,'%03d')])
            end
        end
    end
    ampMap=abs(xMap+1i*yMap);
    file=fileName(ss).name;
    file=[file(1:end-length(fileType)),'mat'];
    cd(resultPathX) %#ok<*MCCD>
    save(file,'xMap')
    cd(resultPathY)
    save(file,'yMap')
%     cd(resultPathCC)
%     save(file,'ccMap')
    cd(resultPathAMP)
    save(file,'ampMap')
    t2=cputime;
    disp(['consuming time',num2str(t2-t1)])
end
    close(h_wait)
    matlabpool close
function pushbutton121_Callback(hObject, eventdata, handles)
uiwait(msgbox('make sure you are in folder of xMap; choose the folder of yMap','Title','modal'));
quiverSize=eval(get(handles.et_quiverSize,'string'));
fileName=handles.fileName;
filePath=handles.filePath;
file=handles.file;
p=length(fileName);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
fileType=handles.fileType;
fileTypeLength=length(fileType);
[file2,filePath2]=uigetfile(file);
cd(filePath2)
fileName2=fileName;
%% binning size
xBinning=eval(get(handles.et_xBinning,'string'));
yBinning=eval(get(handles.et_yBinning,'string'));
%% frame rate unit HZ
frameRate=1000;
%% ROI
xStart=eval(get(handles.et_xBinningStart,'string'));
xEnd=eval(get(handles.et_xBinningEnd,'string'));
yStart=eval(get(handles.et_yBinningStart,'string'));
yEnd=eval(get(handles.et_yBinningEnd,'string'));
%% verify whether ROI fits
if mod(yEnd-yStart+1,yBinning)~=0
    yEnd=floor((yEnd-yStart+1)/yBinning)*yBinning-1+yStart;
end
if mod(xEnd-xStart+1,xBinning)~=0
    xEnd=floor((xEnd-xStart+1)/xBinning)*xBinning-1+xStart;
end
new_xsize=(xEnd-xStart+1)/xBinning; % image size after pixel binning
new_ysize=(yEnd-yStart+1)/yBinning; % image size after pixel binning
%% file name for result
name=['phase_x',num2str(xBinning),'y',num2str(yBinning)];
name=[name,'_',get(handles.et_xy_noise,'string'),'_',num2str(quiverSize(1)),'_',num2str(quiverSize(2))];
result0='phase';
IOSresultPath=get(handles.et_IOSresultPath,'string');
resultPath0=fullfile(IOSresultPath,result0);
cd(IOSresultPath)
if exist(resultPath0)==7
else
    mkdir(result0)
end
cd(resultPath0)
%% subdirectory
result=name;
resultPath=fullfile(resultPath0,result);
cd(resultPath0)
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
[fileTmp,resultPath]=uiputfile('phaseH.mat');
%%
cd(filePath)
img=differentTypeRead(handles.file,fileType);
[m1,n1,o1]=size(img);
handles.newIntensity=inline(get(handles.et_inline,'string'));
newImg=handles.newIntensity(img);
newImg=uint8(newImg);
% figure(90);close(90);figure(90);imshow(newImg);hold on;
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
for ii2=1:p
    %% imread xMap, yMap
    cd(filePath)
    xMap=differentTypeRead(fileName(ii2).name,fileType);
    cd(filePath2)
    yMap=differentTypeRead(fileName2(ii2).name,fileType);
    % phaseH.xMap=xMap;
    % phaseH.yMap=yMap;
    % xMap=xMap2(yStart:yEnd,xStart:xEnd);
    % yMap=yMap2(yStart:yEnd,xStart:xEnd);
    % handles.phaseH_slice.xMap=xMap;
    % handles.phaseH_slice.yMap=yMap;
    %% normalization vectors;
    amplitude=abs(xMap+1i*yMap);
    ampLogic=logical(amplitude);
    xMap(ampLogic)=xMap(ampLogic)./(amplitude(ampLogic));
    yMap(ampLogic)=yMap(ampLogic)./(amplitude(ampLogic));
    xyMap=xMap+1i*yMap;
    %%
    xSect=1:xBinning:(1+new_xsize)*xBinning;
    xSect=xSect+xStart-1;
    xSect(end)=xSect(end)-1;
    ySect=1:yBinning:(1+new_ysize)*yBinning;
    ySect=ySect+yStart-1;
    ySect(end)=ySect(end)-1;
    cx=cell(new_ysize,new_xsize);
    cy=cell(new_ysize,new_xsize);
    phaseMapX=zeros(new_ysize,new_xsize);
    phaseMapY=zeros(new_ysize,new_xsize);
    phaseMap=zeros(new_ysize,new_xsize);
    edgeStep=0;
    for ii=1:new_ysize
        for jj=1:new_xsize
            cx{ii,jj}=[xSect(jj),xSect(jj),xSect(jj+1)-edgeStep,xSect(jj+1)-edgeStep,xSect(jj)].';
            cy{ii,jj}=[ySect(ii),ySect(ii+1)-edgeStep,ySect(ii+1)-edgeStep,ySect(ii),ySect(ii)].';
            phaseMapY(ii,jj)=0.5*(ySect(ii)+ySect(ii+1));
             phaseMapX(ii,jj)=0.5*(xSect(jj)+xSect(jj+1));
     %       plot(cx{ii,jj},cy{ii,jj},'r-','linewidth',1)
    %         plot(phaseMapX(ii,jj),phaseMapY(ii,jj),'r+','linewidth',1)
        end
    end
    %%
    yyPositionTmp=(1:m1).';
    yyPosition=yyPositionTmp*ones(1,n1);
    xxPositionTmp=(1:n1).';
    xxPosition=xxPositionTmp*ones(1,m1);
    xxPosition=xxPosition.';
    mapAll=zeros(m1,n1);
    for ii=1:new_ysize
        for jj=1:new_xsize
            cx2=cx{ii,jj};
            cy2=cy{ii,jj};
            bw2=roipoly(xMap(:,:,1),cx2,cy2);
    %         bwTmp=and(bw,bw2);
            bwTmp=bw2;
            bwXY=and(bwTmp,ampLogic);
            xy=xyMap(bwXY);
            if sum(bwXY(:))>eval(get(handles.et_xy_noise,'string'))
                phaseMap(ii,jj)=sum(xy);
            end
        if abs(phaseMap(ii,jj))~=0
            phaseMap(ii,jj)=phaseMap(ii,jj)/abs(phaseMap(ii,jj));
    %         quiverH=quiver(phaseMapX(ii,jj),phaseMapY(ii,jj),real(phaseMap(ii,jj)),imag(phaseMap(ii,jj)),12);
    %         set(quiverH,'lineWidth',1,'color',[0,1,0])
    %         adjust_quiver_arrowhead_size(quiverH,10);
            mapAll(bw2)=mean(amplitude(bwXY));
        end
        end
    end
    mapAllTmp=mapAll(mapAll>0);
    %%
resizeFactor=2;
if resizeFactor==2
% upperLim=2*mean(mapAllTmp);
    upperLim=1;
    % figure(91);close(91);figure(91);imshow(im2uint8(mat2gray(mapAll,[-upperLim,upperLim])));hold on;
    mapAll=mapAll/2;
    figure(91);close(91);figure(91);imshow(mapAll,[0,upperLim]);hold on;
%         title(['displacement map', ...
%          ' at t = ', num2str(timeCourse(ii2),'%.3f'),' Sec'],'interpreter','none','Fontsize',12)
    map0=jet(256*2);
    map=map0(256:end,:);
    set(gcf,'colormap',map)
    hh=colorbar;
    set(hh,'yTick',linspace(0,upperLim,4))
    % set(hh,'YTick',linspace(0,2,4))
    set(hh,'YTickLabel',{'0','0.1','0.2','0.3'})
    set(hh,'FontSize',15)
else
    upperLim=2;
    % figure(91);close(91);figure(91);imshow(im2uint8(mat2gray(mapAll,[-upperLim,upperLim])));hold on;
%     mapAll=mapAll/2;
    figure(91);close(91);figure(91);imshow(mapAll,[0,upperLim]);hold on;
    figure(91);close(91);figure(91);imshow(mapAll,[0,upperLim]);hold on;
%         title(['displacement map', ...
%          ' at t=', num2str(timeCourse(ii2),'%.3f'),' s'],'interpreter','none','Fontsize',12)
    map0=jet(256*2);
    map=map0(256:end,:);
    set(gcf,'colormap',map)
    hh=colorbar;
    set(hh,'yTick',linspace(0,upperLim,4))
    % set(hh,'YTick',linspace(0,2,4))
    set(hh,'YTickLabel',{'0','0.2','0.4','0.6'})
    set(hh,'FontSize',15)
end
%     %%
%     upperLim=3;
%     % upperLim=2*mean(mapAllTmp);
%     % figure(91);close(91);figure(91);imshow(im2uint8(mat2gray(mapAll,[-upperLim,upperLim])));hold on;
%     figure(91);hold off; imshow(mapAll,[0,upperLim]);hold on;
%     title(['displacement map for image ',fileName(ii2).name(1:end-length(fileType)-1), ...
%          'at', num2str(timeCourse(ii2),'%.3f'),' s'],'interpreter','none')
%     map0=jet(256*2);
%     map=map0(256:end,:);
%     set(gcf,'colormap',map)
%     hh=colorbar;
%     set(hh,'yTick',linspace(0,upperLim,6))
    for ii=1:new_ysize
        for jj=1:new_xsize
            cx2=cx{ii,jj};
            cy2=cy{ii,jj};
            bw2=roipoly(xMap(:,:,1),cx2,cy2);
    %         bwTmp=and(bw,bw2);
            bwTmp=bw2;
            bwXY=and(bwTmp,ampLogic);
            xy=xyMap(bwXY);
            if sum(bwXY(:))>eval(get(handles.et_xy_noise,'string'))
                phaseMap(ii,jj)=sum(xy);
            end
        if abs(phaseMap(ii,jj))~=0
            phaseMap(ii,jj)=phaseMap(ii,jj)/abs(phaseMap(ii,jj));
            theta=angle(phaseMap(ii,jj))+pi;
            phaseMapX(ii,jj)=phaseMapX(ii,jj)+cos(theta)*0.4*xBinning;
            phaseMapY(ii,jj)=phaseMapY(ii,jj)+sin(theta)*0.4*xBinning;            
            quiverH=quiver(phaseMapX(ii,jj),phaseMapY(ii,jj),real(phaseMap(ii,jj)),imag(phaseMap(ii,jj)),quiverSize(1));
            set(quiverH,'lineWidth',2.5,'color',[0,0,1])
            adjust_quiver_arrowhead_size(quiverH,quiverSize(2));
        end
        plot(cx{ii,jj},cy{ii,jj},'r-','linewidth',1)
        end
    end
    rePath_export_figs;
    cd(resultPath)
    colorbar off;
    if get(handles.cb_addT,'value')
        text(580,30,['T=',num2str(timeCourse(ii2),'%2.3f'),' s'],'fontSize',20)
    end
    
    
       set(gca,'color','none')
   export_fig(91,[fileName(ii2).name(1:end-length(fileType)-1),'.png'],'-transparent')
%     saveas(91,fileName(ii2).name(1:end-length(fileType)-1),'tiff')
end
function et_xy_noise_Callback(hObject, eventdata, handles)
% hObject    handle to et_xy_noise (see GCBO)
% et_xy_noise as text
%        str2double(get(hObject,'String')) returns contents of et_xy_noise as a double
function et_xy_noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_xy_noise (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_quiverSize_Callback(hObject, eventdata, handles)
% hObject    handle to et_quiverSize (see GCBO)
% et_quiverSize as text
%        str2double(get(hObject,'String')) returns contents of et_quiverSize as a double
function et_quiverSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_quiverSize (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% --------------------------------------------------------------------
function mp_roiZero_Callback(hObject, eventdata, handles)
% hObject    handle to mp_roiZero (see GCBO)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
result='ROIsetzero';
resultFolder=fullfile(handles.filePath,result);

if exist(resultFolder)==7
else
    cd(filePath)
    mkdir(result)
end
cd(resultFolder)
[fileTmp,resultFolder]=uiputfile(file);
prompt={'background:'};
   name='backgournd';
   numlines=1;
   defaultanswer={'0'};
    options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
   answer=inputdlg(prompt,name,numlines,defaultanswer,options);
   ratio=eval(answer{1});
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
rawData=zeros(p,sum(ROI_N(:))+1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
h_wait=waitbar(0,'wait');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    cd(filePath)
    img=differentTypeRead(file,fileType);
    vv=1;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
            vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            img(bw)=ratio;
%             rawData(kk,vv)=mean(ROItmp(:));
        end
    end
    cd(resultFolder)
    differentTypeWrite(img,file,fileType,handles.imgDepth);
end
close(h_wait)
rawData(:,2:end)= ...
    IOS_time_gui_filter(rawData(:,2:end),handles,'vertical');
% rawData=log(rawData);
handles.rawData_meanOfROI=rawData;
guidata(hObject,handles)
function et_movie_fps_Callback(hObject, eventdata, handles)
function et_movie_fps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in cb_movie_curve.
function cb_movie_curve_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    if isfield(handles,'movie_curve')
    else
        handles.movie_curve=1;
        file=handles.file;
        filePath=handles.filePath;
        fileType=handles.fileType;
        cd(filePath)
        img=differentTypeRead(file,fileType);
        [mm,nn,nn3]=size(img);
        set(handles.et_movie_range,'string',['[ 1:', num2str(nn),' ]'])
    end
end
p=handles.p;
stringT=['1:1:',num2str(p)];
set(handles.et_movie_p,'string',stringT)
guidata(hObject,handles)
function et_movie_range_Callback(hObject, eventdata, handles)
% hObject    handle to et_movie_range (see GCBO)
% et_movie_range as text
%        str2double(get(hObject,'String')) returns contents of et_movie_range as a double
function et_movie_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_movie_range (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pb_movie_run.
function pb_movie_run_Callback(hObject, eventdata, handles)
% hObject    handle to pb_movie_run (see GCBO)
function et_movie_y_Callback(hObject, eventdata, handles)
% hObject    handle to et_movie_y (see GCBO)
% et_movie_y as text
%        str2double(get(hObject,'String')) returns contents of et_movie_y as a double
function et_movie_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_movie_y (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton124_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton124 (see GCBO)
function pushbutton125_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton125 (see GCBO)
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
%%
filePathTmp=get(handles.et_IOSresultPath,'string');
resultName1='movie';
resultPath1=fullfile(filePathTmp,resultName1);
cd(filePathTmp)
if exist(resultPath1)==7
else
    mkdir(resultName1);
end
cd(filePath)
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
fileName=handles.fileName;
p=length(fileName);
newIntensity=inline(get(handles.et_inline,'string'));
contents = cellstr(get(handles.pm_colormap_raw,'String'));
colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
colorSelected=eval(colorSelectedTmp{1});
%%
bscan=handles.bscan;
cd(filePath)
img=differentTypeRead(file,fileType);
[mm,nn,n3]=size(img);
bscan0=uint8(255*ind2rgb(round(newIntensity(bscan)),colorSelected));
pp2=size(bscan0,1);
%%
singleFrame=get(handles.cb_movie_single,'value');
stringT=eval(get(handles.et_movie_p,'string'));
crossY=eval(get(handles.et_movie_y,'string'));
%%
hh=waitbar(0,'please wait ...');
 cd(filePath)
 rePath_export_figs;






 if get(handles.cb_movie_curve,'value')
     figure(84);
    imgAll=255*ones(mm+p+10,nn,3,'uint8');
     aa=eval(get(handles.et_movie_range,'string'));
     curveB=mean(bscan0(:,aa),2);
     x=mm+11:mm+p+10;x2=x;
     bas_amp=eval(get(handles.et_movie_bas_amp,'string'));
     for ii=1:p
         figure(84);hold off;
                    waitbar(ii/p,hh, ...
                     [num2str(100*ii/p,'%03.1f'),'% completed']);
                 cd(filePath)
        file=fileName(ii).name;
        img=differentTypeRead(file,fileType);
        %% adjust intensity
      newImg=handles.newIntensity(img);
        if get(handles.cb_eval,'value')
            x=newImg;
            eval(get(handles.et_eval,'string')) ;
            newImg=x;
        end
        
        img=newImg;
        %% adjust itensity done
%         img2=newIntensity(img);
        img2=uint8(255*ind2rgb(round(newIntensity(img)),colorSelected));
        img2(crossY,1:1,1)=255;
        img2(crossY,1:1,2)=0;
        img2(crossY,1:1,3)=0;
%         img2(crossY(1),1:10,1)=255;
%         img2(crossY(1),1:10,2)=0;
%         img2(crossY(1),1:10,3)=0;
%         img2(crossY(end),1:10,1)=255;
%         img2(crossY(end),1:10,2)=0;
%         img2(crossY(end),1:10,3)=0;
        imgAll(1:mm,1:nn,:)=img2;
        bscan1=bscan0;
        bscan1(stringT(ii),:,1)=255;
        bscan1(stringT(ii),:,2)=0;
        bscan1(stringT(ii),:,3)=0;
        imgAll(mm+10+1:mm+10+pp2,:,:)=bscan1;
        figure(84); imshow(uint8(imgAll));
        hold on;
        plot(bas_amp(1)+bas_amp(2)*curveB,x2,'r')
        if singleFrame
                cd(resultPath1)
                saveas(gcf,file(1:end-length(fileType)-1),'tiff')
%                    set(gca,'color','none')
%                export_fig(gcf,[file(1:end-length(fileType)-1),'.png'],'-transparent')
%         imwrite(imgAll,[file(1:end-length(fileType)),'png']);
        end
        %A(:,i)=getframe(fig1,winsize);
        %http://www.math.canterbury.ac.nz/~c.scarrott/MATLAB_Movies/movies.html
%         [image1Ind,colorMap_1]=rgb2ind(imgAll,256);
%         M_color(ii)=im2frame(image1Ind,colorMap_1);
            M_color(ii)=getframe(gcf);
    end
 else
imgAll=255*ones(mm+p+10,nn,3,'uint8');
     for ii=1:p
                    waitbar(ii/p,hh, ...
                     [num2str(100*ii/p,'%03.1f'),'% completed']);
                 cd(filePath)
        file=fileName(ii).name;
        img=differentTypeRead(file,fileType);
        %% adjust intensity
      newImg=handles.newIntensity(img);
        if get(handles.cb_eval,'value')
            x=newImg;
            eval(get(handles.et_eval,'string')) ;
            newImg=x;
        end
        img=newImg;
        %% adjust itensity done
%         img2=newIntensity(img);
        img2=uint8(255*ind2rgb(round(newIntensity(img)),colorSelected));
        img2(crossY,1:1,1)=255;
        img2(crossY,1:1,2)=0;
        img2(crossY,1:1,3)=0;
        img2(crossY(1),1:10,1)=255;
        img2(crossY(1),1:10,2)=0;
        img2(crossY(1),1:10,3)=0;
        img2(crossY(end),1:10,1)=255;
        img2(crossY(end),1:10,2)=0;
        img2(crossY(end),1:10,3)=0;
        imgAll(1:mm,1:nn,:)=img2;
        bscan1=bscan0;
        bscan1(stringT(ii),:,1)=255;
        bscan1(stringT(ii),:,2)=0;
        bscan1(stringT(ii),:,3)=0;
        imgAll(mm+10+1:mm+10+pp2,:,:)=bscan1;
        if singleFrame
                cd(resultPath1)
        imwrite(imgAll,[file(1:end-length(fileType)),'png']);
        end
%         [image1Ind,colorMap_1]=rgb2ind(imgAll,256);
%         M_color(ii)=im2frame(image1Ind,colorMap_1);
            M_color(ii)=im2frame(imgAll,colorSelected);
    end
 end
close(hh)
    cd(resultPath1)
     imwrite(uint8(bscan),[file(1:end-length(fileType)),'tif'])
[movieName,resultPath1]=uiputfile('movie.avi');
movieName=movieName(1:end-4);
cd(resultPath1)
movieFPS=eval(get(handles.et_movie_fps,'string'));
hh=waitbar(.50,'saveing raw. please wait ...');
 movie2avi(M_color,fullfile(resultPath1,movieName),'compression','none','fps',movieFPS);
 close(hh)
% movie2avi(M_color,fullfile(resultPath1,movieName),'compression','MSVC','quality',85,'fps',movieFPS);
%     'compression','MSVC'
% movie2avi(M_color,fullfile(resultPath1,movieName),'compression','none','fps',movieFPS);
    % movie2avi(M2,fullfile(filePath,'movie2'),'compression','none','fps',5);
%     implay(M_color,movieFPS)
guidata(hObject,handles)
function pushbutton126_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton126 (see GCBO)
function pushbutton127_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton127 (see GCBO)
file=handles.file;
fileType=handles.fileType;
filePath=handles.filePath;
%%
cd(filePath)
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
fileName=handles.fileName;
p=length(fileName);
newIntensity=inline(get(handles.et_inline,'string'));
contents = cellstr(get(handles.pm_colormap_raw,'String'));
colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
colorSelected=eval(colorSelectedTmp{1});
%%
bscanP=eval(get(handles.et_movie_y,'string'));
cd(filePath)
img=differentTypeRead(file,fileType);
[mm,nn,n3]=size(img);
bscan=zeros(p,nn);
hh=waitbar(0,'please wait ...');
 for ii=1:p
                waitbar(ii/p,hh, ...
                 [num2str(100*ii/p,'%03.1f'),'% completed'],'get bscan');
    file=fileName(ii).name;
    img=differentTypeRead(file,fileType);
    bscan(ii,:)=mean(img(bscanP,:,1));
 end
close(hh)
figure(71); hold off; imshow(uint8(newIntensity(bscan)),colorSelected);
hold on;
aa=eval(get(handles.et_movie_range,'string'));
curveB=mean(bscan(:,aa),2);
x=1:p;
 bas_amp=eval(get(handles.et_movie_bas_amp,'string'));
%  plot(x,bas_amp(1)+bas_amp(2)*curveB,'r')
plot(bas_amp(1)+bas_amp(2)*curveB,x,'r')
handles.bscan=bscan;
handles.bscan2=bscan;
% bscanP=eval(get(handles.et_movie_y,'string'));
cd(filePath)
result='bscan';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
resultName=['bscan',num2str(bscanP(1),'%03d'),'_',num2str(bscanP(end),'%03d')];
[resultName,resultPath]=uiputfile(resultName,'save bscan');
resultName1=[resultName,'_uint8.tif'];
resultName2=[resultName,'_uint16.tif'];
cd(resultPath)
imwrite(uint8(bscan),resultName1,'tif')
imwrite(uint16(bscan),resultName2,'tif')
guidata(hObject,handles)
function cb_movie_single_Callback(hObject, eventdata, handles)
function et_movie_compress_Callback(hObject, eventdata, handles)
function et_movie_compress_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text118_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.*','select movie please');
cd(filePath)
aviObject=VideoReader(file);
p=aviObject.NumberOfFrames;
dotNO=strfind(file,'.');
dotNO=dotNO(end);
fileType=file(dotNO+1:end);
fileTmp=file(1:dotNO-1);
%% see the first frame
image1=read(aviObject,1);
figure(3); imshow(image1,[]);
[image1Ind,colorMap_1]=rgb2ind(image1,256);
[mm,nn]=size(image1(:,:,1));
x1=1;y1=1;
x2=nn;y2=mm;
% x1=288; y1=1;
% x2=930; y2=835;
% x1=120; y1=189;
% x2=1082;y2=800;
% x1=231; x2=894;y1=1;y2=mm;
h_wait=waitbar(0,'please wait');
for ii=1:p
    waitbar(ii/p,h_wait,[num2str(ii*100/p,'%03.1f'),'completed']);
    image1=read(aviObject,ii);
    image1Ind=rgb2ind(image1,colorMap_1);
    X=image1Ind(y1:y2,x1:x2);
    MAP=colorMap_1;
     M(ii)=im2frame(X,MAP);
%      imwrite(image1,[num2str(ii,'%03d'),'.tif'],'tif');
end
close(h_wait)
% implay(M)
% movie2avi(M,fullfile(filePath,'movieCompression'),'compression','none','fps',16);
% movie2avi(M,fullfile(filePath,'movieCompressionMSVCQuality75'),'compression','MSVC','fps',16);
%movie2avi(M,fullfile(filePath,'movieCompressionMSVCQuality85'),'compression','MSVC','quality',85,'fps',16);
movie2avi(M,fullfile(filePath,[fileTmp,'_cprs', get(handles.et_movie_compress,'string')]), ...
    'compression','MSVC', ...
    'quality',eval(get(handles.et_movie_compress,'string')), ...
    'fps',aviObject.FrameRate);
function et_movie_bas_amp_Callback(hObject, eventdata, handles)
function et_movie_bas_amp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function mp_removeBackground_Callback(hObject, eventdata, handles)
cd(handles.filePath)
%% imread 1
file1=handles.file;
filePath1=handles.filePath;
% [file1,filePath1]=uigetfile('*.*','select 1st series of images');
% cd(filePath1);
indexNO=strfind(file1,'.');
fileType1=file1(indexNO(end)+1:end);
% fileName1=dir(['*.',fileType1]);
fileName1=handles.fileName;
p=length(fileName1);cd(filePath1)
img=differentTypeRead(file1,fileType1);
% img=imread(file1);
[m1,n1,q1]=size(img);
img3=zeros(size(img));
%% imread 2
[file2,filePath2]=uigetfile('*.*','select 2n series of images');
cd(filePath2);
indexNO=strfind(file2,'.');
fileType2=file2(indexNO(end)+1:end);
fileName2=dir(['*.',fileType2]);
img=differentTypeRead(file2,fileType2);
[m2,n2,q2]=size(img);
cd(filePath1)
result='IOS_pure';
resultFolder=fullfile(filePath1,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
[result,resultPath]=uiputfile(file1,'choose a folder to save');
%%
 prompt={'background:'};
   name='set background (only 1st useful for gray; all for RGB)';
   numlines=1;
   defaultanswer={'[128,256,128]'};
    options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
   answer=inputdlg(prompt,name,numlines,defaultanswer,options);
   ratio=eval(answer{1});
 %  ratio=ratio/sum(ratio);
%%
h_wait=waitbar(0,'wait');
    for ii=1:p
            waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
        cd(filePath1)
        file1=fileName1(ii).name;
        img1=differentTypeRead(file1,fileType1);
        cd(filePath2)
        file2=[file1(1:end-length(fileType1)),fileType2];
        img2=differentTypeRead(file2,fileType2);
        if q1==3 && q2==1
                     img3(:,:,1)=img1(:,:,1).*img2/1;
        img3(:,:,2)=img1(:,:,2).*img2/1;
        img3(:,:,3)=img1(:,:,3).*img2/1;
%          img3(:,:,1)=img1(:,:,1).*img2/max(img2(:));
%         img3(:,:,2)=img1(:,:,2).*img2/max(img2(:));
%         img3(:,:,3)=img1(:,:,3).*img2/max(img2(:));
        imgFlag=~logical(img2);
        imgTmp=img3(:,:,1);
        imgTmp(imgFlag)=ratio(1);
        img3(:,:,1)=imgTmp;
        %
        imgTmp=img3(:,:,2);
        imgTmp(imgFlag)=ratio(2);
        img3(:,:,2)=imgTmp;
        imgTmp=img3(:,:,3);
        imgTmp(imgFlag)=ratio(3);
        img3(:,:,3)=imgTmp;
        elseif q1==1 && q2==1
            img3=img1.*img2;
            img3(~logical(img2))=ratio(1);
        end
        cd(resultPath)
        differentTypeWrite(img3,file1, handles.fileType,handles.imgDepth);
%         imwrite(img1,file1);
    end
close(h_wait)
function t_positions_ButtonDownFcn(hObject, eventdata, handles)
function et_movie_p_Callback(hObject, eventdata, handles)
function et_movie_p_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton128_Callback(hObject, eventdata, handles)
%%
 prompt={'strech b scan:'};
   name='Input for stretch function';
   numlines=1;
   if isfield(handles,'bScanStretch')
       defaultanswer=handles.bScanStretch;
   else
       defaultanswer={'im2uint8(mat2gray(x))'};
   end
    options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
   answer=inputdlg(prompt,name,numlines,defaultanswer,options);
   newIntensity=inline(answer{1});
   handles.bScanStretch=answer(1);
   figure(72); subplot(2,1,1); imshow(uint8(handles.bscan2))
   handles.bscan=newIntensity(handles.bscan2);
   subplot(2,1,2); imshow(uint8(handles.bscan))
   guidata(hObject,handles)
   filePath=handles.filePath;
   bscanP=eval(get(handles.et_movie_y,'string'));
cd(filePath)
result='bscan';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(result)
end
cd(resultPath)
resultName=['bscan',num2str(bscanP(1),'%03d'),'_',num2str(bscanP(end),'%03d'),'_stretch'];
[resultName,resultPath]=uiputfile(resultName,'save bscan');
resultName1=[resultName,'_uint8.tif'];
% resultName2=[resultName,'_uint16.tif'];
cd(resultPath)
imwrite(uint8(handles.bscan),resultName1,'tif')
% imwrite(uint16(bscan),resultName2,'tif')
guidata(hObject,handles)
function cb_vertical_Callback(hObject, eventdata, handles)
function cb_continuous_Callback(hObject, eventdata, handles)
function rePath_export_figs
pathWhich=which('IOS_Software');
dotNO=strfind(pathWhich,filesep);
dotNO=dotNO(end);
pathWhich=pathWhich(1:dotNO-1);
filePath=pathWhich;
file='export_figs';
filePath=fullfile(filePath,file);
path(filePath,path);
function popupmenu18_Callback(hObject, eventdata, handles)
rePath_export_figs;
method=get(handles.popupmenu18,'value');
fig=eval(get(handles.et_exportFig,'string'));
 figure(fig);
if method==1
    [result,resultPath]=uiputfile('saveFigure_transparent.png');
    cd(resultPath)
   set(gca,'color','none')
   export_fig(fig,result,'-transparent')
elseif method==2
    [result,resultPath]=uiputfile('saveFigure_transparent_native.png');
    cd(resultPath)
   set(gca,'color','none')
   export_fig(fig,result,'-transparent','-native')
elseif method==3
    [result,resultPath]=uiputfile('saveFigure_transparent_painter.png');
    cd(resultPath)
   set(gca,'color','none')
   export_fig(fig,result,'-transparent','-painters')
elseif method==4
    [result,resultPath]=uiputfile('saveFigure_transparent_painter.pdf');
    cd(resultPath)
%    set(gca,'color','none')
   export_fig(fig,result,'-q101')
elseif method==5
     open('export_fig - Oliver Woodford.htm');
end
function popupmenu18_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function et_exportFig_Callback(hObject, eventdata, handles)
function et_exportFig_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton129_Callback(hObject, eventdata, handles)
if isfield(handles,'comString')
    comString=handles.comString;
else
    comString='plot(0+1*y,x,plotGroup{ii},''lineWidth'',1)';
end
% handles.comString
 prompt={'combine:'};
   name='Input base and amplitude';
   numlines=1;
   defaultanswer={comString};
    options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
   answer=inputdlg(prompt,name,numlines,defaultanswer,options);
   bas_amp=answer{1};
   handles.comString=bas_amp;
file=handles.file;
filePath=handles.filePath;
fileType=handles.fileType;
cd(filePath)
img=differentTypeRead(file,fileType);
img=img(:,:,1);
if get(handles.cb_vertical,'value')==0
    img=img.';
end
[mm,nn]=size(img(:,:,1));
eval(get(handles.et_verticalLines,'string'));
p=length(yr);
Y=zeros(mm,p);
X=Y;
fig1=31;
fig2=32;
figure(fig1);close(fig1);figure(fig1);hold on;
figure(fig2);close(fig2);figure(fig2);hold on;
plotGroup={'r','g','b','r-*','g-*','b-*','r--','g--','b--'};
titleName='';
newImg=getCurrentImg(handles);
figure(fig2);imshow(uint8(newImg));hold on;
for ii=1:p
    yRange=yr{ii};
    y=mean(img(:,yRange),2);
    y=IOS_time_gui_filter(y,handles,'vertical');
    x=[1:mm].';
    figure(fig1);plot(x,y,plotGroup{ii});grid on;
    Y(:,ii)=y;
    X(:,ii)=x;
    titleName=[titleName,num2str(yRange(1)),' to ',num2str(yRange(end)),';'];
    figure(fig2);eval(bas_amp)
end
title(titleName)
handles.vertical.Y=Y;
handles.vertical.X=X;
handles.vertical.Fig1=fig1;
handles.vertical.Fig2=fig2;
    guidata(hObject,handles)
function newImg=getCurrentImg(handles)
contents = cellstr(get(handles.pm_colormap_raw,'String'));
colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
colorSelected=eval(colorSelectedTmp{1});
cd(handles.filePath)
handles.newIntensity=inline(get(handles.et_inline,'string'));
fileType=handles.fileType;
file=handles.file;
img=differentTypeRead(fullfile(handles.filePath,file),fileType);
 newImg=handles.newIntensity(img);
%% further processing
if get(handles.cb_eval,'value')
    x=newImg;
    eval(get(handles.et_eval,'string')) ;
    newImg=x;
end
function cb_uint16_Callback(hObject, eventdata, handles)
function cb_IOS_continuous_Callback(hObject, eventdata, handles)
function me_fft_Callback(hObject, eventdata, handles)
function Untitled_13_Callback(hObject, eventdata, handles)
cd(handles.filePath)
handles.newIntensity=inline(get(handles.et_inline,'string'));
fileType=handles.fileType;
file=handles.file;
x=differentTypeRead(file,fileType);
d3=size(x,3);
fileName=handles.fileName;
p=length(fileName);
filePath=handles.filePath;
result1='fft2';
resultPath1=fullfile(filePath,result1);
cd(filePath)
if exist(resultPath1)==7
else
    mkdir(result1)
end
cd(resultPath1)
[file,resultPath1]=uiputfile(file);
h_wait=waitbar(0,'wait');
for ii=1:p
    waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed']);
%     cd(filePath)
    file1=fileName(ii).name;
    x=differentTypeRead(fullfile(filePath,file1),fileType);
    if max(x(:))>0
    y=fftshift(fft2(x));
    else
        y=0*x;
    end
    file1=[file1(1:end-length(fileType)),'mat'];
    save(fullfile(resultPath1,file1),'y')     
end 
  close(h_wait)


% --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
function me_center_Callback(hObject, eventdata, handles)
if isfield(handles,'centerR')
    centerR=handles.centerR;
else
    centerR='[104,104,50,0.6]';
end
prompt={'input x,y,r,pixelSize'};
name='Input x, y, r';
numlines=1;
defaultanswer={centerR};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer=inputdlg(prompt,name,numlines,defaultanswer,options);
handles.centerR=answer{1};
ratio=eval(answer{1});
x0=ratio(1);
y0=ratio(2);
r=ratio(3);
pixelSize=ratio(4);
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
img=differentTypeRead(fullfile(filePath,file),fileType);
img=img(:,:,1);
[m,n]=size(img);
a=[x0,n-x0+1,y0,n-y0+1,r];
r=round(min(a));
y=zeros(r,2);
y(:,1)=[1:r].';
y(:,1)=y(:,1)*(1/pixelSize/m*10);
y(1,2)=img(y0,x0);

y1=1:m;
x1=1:n;
[xx1,yy1]=meshgrid(x1,y1);
h_wait=waitbar(0,'wait');
for ii=r-1:-1:1
    waitbar(ii/r,h_wait,[num2str(100*ii/r,'%04.1f'),'%left']);
    rMax=ii;
    thetaTmp=linspace(0,360,round(2*pi*rMax));
    thetaY=y0+sin(thetaTmp)*rMax;
    thetaX=x0+cos(thetaTmp)*rMax;
    vq=griddata(xx1,yy1,double(img),thetaX,thetaY,'nearest');
    y(ii+1,2)=mean(vq(:));
end
close(h_wait);
figure(6);plot(y(:,1),y(:,2),'r-*')
xlabel('cycles/10\mum');
grid on;
handles.centerY=y;
guidata(hObject,handles)

function Untitled_17_Callback(hObject, eventdata, handles)
cd(handles.filePath)
centerY=handles.centerY;
[file,filePath]=uiputfile('centerToSurround.txt');
cd(filePath)
save(file,'centerY','-ASCII');
saveas(6,'centerY','tiff');

function cb_dF_Callback(hObject, eventdata, handles)

function cb_map_raw_Callback(hObject, eventdata, handles)

function cb_addT_Callback(hObject, eventdata, handles)

function pb_filterTime_Callback(hObject, eventdata, handles)


file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
imgDepth=handles.imgDepth;
[m,n,ntrash]=size(img);

%%
nameTmp=handles.filterMethod.name;
rawData=zeros(p,1);
if nameTmp==4
    h_wait=waitbar(0,'wait');
    for kk=1:p
        %h_wait=waitbar(0,'wait');
        waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed: getting parameters from detrending G1R1']);
        file=fileName(kk).name;
        img=differentTypeRead(fullfile(filePath,file),fileType);
%         vv=1;
        for ii=1:1
            for jj=1:1
%                 vv=vv+1;
                cx=handles.pointPst{ii}{jj}(:,1);
                cy=handles.pointPst{ii}{jj}(:,2);
                bw=roipoly(img(:,:,1),cx,cy);
                ROItmp=img(bw);
                rawData(kk)=mean(ROItmp(:));
            end
        end
    end
    close(h_wait)
    detrend_order=handles.filterMethod.detrend_order;
    detrend_section=handles.filterMethod.detrend_section;
    detrend_base=handles.filterMethod.detrend_base;
    if strcmp(detrend_base,'mean')
        t1=detrendnonlin(rawData,detrend_order,detrend_section);
    else
        t1=detrendnonlin(rawData,detrend_order,detrend_section,eval(detrend_base));
    end
    
    background=(rawData-t1);
    
    h_wait=waitbar(0,'wait');
    result=['detrend_order',num2str(detrend_order)];
    resultPath=fullfile(filePath,result);
    if exist(resultPath)==7
        
    else
        cd(filePath)
        mkdir(result)
    end
     for kk=1:p
        waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed detrending images']);
        file=fileName(kk).name;
        img=differentTypeRead(fullfile(filePath,file),fileType);
        img=img-background(kk);
        resultFile=[file(1:end-length(fileType)),'mat'];
        save(fullfile(resultPath,resultFile),'img');
             
    end
    close(h_wait)
    
else

    if strcmp(fileType,'mat')
        imgStore=zeros(p,m*n,'single');
    else
        imgStore=zeros(p,m*n,imgDepth);

    end
    %% get img
    h_wait=waitbar(0,'wait');
    for ii=1:p
        %h_wait=waitbar(0,'wait');
        waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed for reading images']);    
        file=fileName(ii).name;
        if strcmp(fileType,'mat')
            img=differentTypeRead(fullfile(filePath,file),fileType);
            imgStore(ii,:)=img(:).';
        else
            img2=imread(fullfile(filePath,file));
            imgStore(ii,:)=img2(:).';

        end

    end
    %% filter temporally
    for ii=1:m*n
        waitbar(ii/m/n,h_wait,[num2str(100*ii/m/n,'%04.1f'),'%completed for filtering images']);    
        img1Point=imgStore(:,ii);
            img1Point2=IOS_time_gui_filter(single(img1Point),handles,'vertical');
        if strcmp(fileType,'mat')

            imgStore(:,ii)=img1Point2;
        else
            imgStore(:,ii)=eval([imgDepth,'(img1Point2)']);
        end
    end
    %% save img
    result='fltrTmpry';
    resultPath=fullfile(filePath,result);
    pathDe=cd;
    cd(filePath)
    if exist(resultPath)==7
    else
        mkdir(result)
    end
    cd(resultPath)
    [file,resultPath]=uiputfile(file);
    cd(pathDe)
    for ii=1:p
        waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed for saving images']);    
        file=fileName(ii).name;
        img3=reshape(imgStore(ii,:).',m,n);
        differentTypeWrite(img3,fullfile(resultPath,file),fileType,imgDepth);

    end
    close(h_wait)    
end


function me_ROIauto_Callback(hObject, eventdata, handles)




function differentTypeShow(img,fig,handles)
contents = cellstr(get(handles.pm_colormap_raw,'String'));
colorSelectedTmp=contents(get(handles.pm_colormap_raw,'value'));
if get(handles.cb_uint16,'value')
colorSelectedTmp{1}=[colorSelectedTmp{1}(1:end-5),'(2^16)'];
end
colorSelected=eval(colorSelectedTmp{1});
handles.newIntensity=inline(get(handles.et_inline,'string'));
 newImg=handles.newIntensity(img);
%% further processing
if get(handles.cb_eval,'value')
    x=newImg;
    eval(get(handles.et_eval,'string')) ;
    newImg=x;
end
figure(fig);
if get(handles.cb_uint16,'value')
    newImgColor=ind2rgb(uint16(newImg),colorSelected);
else
    newImgColor=ind2rgb(uint8(newImg),colorSelected);
end  
imshow(newImgColor)


function me_ROIseperate_Callback(hObject, eventdata, handles)
prompt = {'select ROI with interval:','select searching time zone'};
dlg_title = 'Input for selecting ROI';
num_lines = 1;
if isfield(handles,'ROIintDef')
    def=handles.ROIintDef;
else
    def = {'15',['1:',num2str(handles.p)]};
end


options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
handles.ROIintDef=answer;
interval=eval(answer{1});
searZone=eval(answer{2});
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
fig=15;
differentTypeShow(img,fig,handles);
figure(fig);
hold on;
[mm,nn,nntrash]=size(img);
x=1:interval:nn;
px=length(x)-1;
sect=cell(px,1);
for ii=1:px
    sect{ii}=x(ii):x(ii+1)-1;
    plot([x(ii),x(ii)],[1,mm],'r')
end
ii=ii+1;
plot([x(ii),x(ii)],[1,mm],'r')
hold off;
title(['interval=',num2str(interval)])
data=zeros(p,px);
h_wait=waitbar(0,'please wait');

for ii=1:p
    waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed for saving images']);  
    file=fileName(ii).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    for jj=1:px
        imgTmp=img(1:mm,sect{jj});
        data(ii,jj)=mean(imgTmp(:));
    end
    
end
close(h_wait)
data(:,1:end)= ...
    IOS_time_gui_filter(data(:,1:end),handles,'vertical');
% dataMean=mean(data);
data2=data;
pre=eval(get(handles.et_pre2,'string'));
for ii=1:px
    data2(:,ii)=(data2(:,ii)-mean(data2(1:pre,ii)))/mean(data2(1:pre,ii));
end

%% get statistics
dataMean=mean(data2(1:pre,:));
dataStd=std(data2(1:pre,:));
[dataMin,Ind]=min(data2(searZone,:));
Ind=Ind+searZone(1)-1;
minAll=0;
figure(16);hold off;grid on;
for ii=1:px
    dataPlot=data2(:,ii);
    dataPlot=dataPlot+minAll;
    if ii==2
        hold on;
    end
        plot(1:p,dataPlot,Ind(ii),dataPlot(Ind(ii)),'o')
    minAll=minAll+dataMin(ii)*1.5;
    
    
end
title('IOS curve')
figure(17); hold off; plot(1:length(dataMin),abs(dataMin),'-o');
grid on;
title('IOS')
guidata(hObject,handles)

%% get statistics
data3=data;
for ii=1:px
    data3(:,ii)=(data3(:,ii)-mean(data3(1:pre,ii)));
end
[dataMin3,Ind]=min(data3(searZone,:));
figure(18); hold off; plot(1:length(dataMin3),abs(dataMin3),'-o');
grid on;
title('abslute')
result='lobsterAngle';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end
saveas(18,fullfile(resultPath,'absolute.tif'),'tiff')
saveas(17,fullfile(resultPath,'IOS.tif'),'tiff')
saveas(16,fullfile(resultPath,'IOSraw.tif'),'tiff')
saveas(15,fullfile(resultPath,'Imgraw.tif'),'tiff')
save(fullfile(resultPath,'absolute.txt'),'dataMin3','-ASCII')
save(fullfile(resultPath,'IOS.txt'),'dataMin','-ASCII')
save(fullfile(resultPath,'IOS_raw.txt'),'data','-ASCII')

function me_ROIautoContinue_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
[mm,nn,nntrash]=size(img);

prompt = {'select ROI with interval:','select searching time zone','center point','how many degree/pixel'};
dlg_title = 'Input for selecting ROI';
num_lines = 1;
if isfield(handles,'ROIintDefc')
    defc=handles.ROIintDefc;
else
    defc = {'15',['1:',num2str(handles.p)],num2str(round(nn/2)),'0.1767'};
end


options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer = inputdlg(prompt,dlg_title,num_lines,defc,options);
handles.ROIintDefc=answer;
interval=eval(answer{1});
searZone=eval(answer{2});
centerP=eval(answer{3});
unit=eval(answer{4});

%% reset sect continuously
px=nn-interval+1;
sect=cell(px,1);
for ii=1:px
    sect{ii}=(ii-1)+1:(ii-1)+interval;
end

data=zeros(p,px);
h_wait=waitbar(0,'please wait');

for ii=1:p
    waitbar(ii/p,h_wait,[num2str(100*ii/p,'%04.1f'),'%completed for saving images']);  
    file=fileName(ii).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    for jj=1:px
        imgTmp=img(1:mm,sect{jj});
        data(ii,jj)=mean(imgTmp(:));
    end
    
end
close(h_wait)
data(:,1:end)= ...
    IOS_time_gui_filter(data(:,1:end),handles,'vertical');
handles.ROIcon.data=data;
handles.ROIcon.answer=answer;
guidata(hObject,handles)
continousROIsShow(handles);
guidata(hObject,handles)
% dataMean=mean(data);

%%
function continousROIsShow(handles)
%%
ROIcon=handles.ROIcon;
data=ROIcon.data;
answer=ROIcon.answer;

fig=15;
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
[mm,nn,nntrash]=size(img);

differentTypeShow(img,fig,handles);
figure(fig);
hold on;
interval=eval(answer{1});
searZone=eval(answer{2});
centerP=eval(answer{3});
unit=eval(answer{4});
px=nn-interval+1;




data2=data;
pre=eval(get(handles.et_pre2,'string'));
for ii=1:px
    data2(:,ii)=(data2(:,ii)-mean(data2(1:pre,ii)))/mean(data2(1:pre,ii));
end

%% get statistics
centerP=eval(answer{3});
centerP=centerP-floor(interval/2);
centerP1=centerP:-interval:1;
centerP1=centerP1(end:-1:1);
centerP2=centerP:interval:px;
centerP3=[centerP1(1:end-1),centerP2];


%%
px2=length(centerP3);
for ii=1:px2

    plot([centerP3(ii),centerP3(ii)],[1,mm],'r')
end

plot([centerP3(ii)+interval-1,centerP3(ii)+interval-1],[1,mm],'r')
hold off;
title(['interval=',num2str(interval)])

%%
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';

dataMean=mean(data2(1:pre,:));
dataStd=std(data2(1:pre,:));
[dataMin,Ind]=min(data2(searZone,:));
% dataMin=mean(data2(450:600,:));
Ind=Ind+searZone(1)-1;
minAll=0;

%% dealign with half peak
figure(19); hold off; 
x=1:px;x=(x-centerP)*unit;
halfMin=dataMin/2;
halfTime=zeros(1,px);
halfTime2=zeros(1,px);
maxValueFitting=zeros(px,1);
xTime=1:p;
for ii=1:px
    halfData2=data2(searZone,ii);
    if ii==1
        [timeToHalf,maxValue]=findHalf_expo(timeCourse,data2(:,ii),searZone(1),searZone(1)+200,'fig');
    else
       [timeToHalf,maxValue]=findHalf_expo(timeCourse,data2(:,ii),searZone(1),searZone(1)+200); 
    end    
    maxValueFitting(ii)=maxValue;
%     if ii==62
%         trash=0;
%     end
    ps=length(searZone);
    r1=1:ps;
    r2=1:0.001:ps;
    halfData3=interp1(r1,halfData2,r2,'linear');

    
    
%     halfData=data2(searZone,ii);
%     yi = interp1(x,y,xi);
    halfData=abs(halfData3-halfMin(ii));
    IndTmp=find(halfData<0.0001);
    if isempty(IndTmp)
        IndTmp=5;
    end
    IndTmp=IndTmp(1);
    %%
%     x1=IndTmp;
%     
%     if halfData2(x1)-halfMin(ii)>=0
%         x2=x1+1;
%     else
%         x2=x1-1;
%     end
%     y1=halfData2(x1);
%     y2=halfData2(x2);
%     y0=halfMin(ii);
%     x0=(y0-y1)*(x2-x1)/(y2-y1)+x1;
 
    halfTime(ii)=r2(IndTmp)+searZone(1)-1;
    halfTime(ii)=interp1(xTime,timeCourse,halfTime(ii),'linear');
    halfTime2(ii)=timeToHalf;
end
figure(19);
 subplot(2,1,1);
 plot(x,halfTime,'-*');xlabel('Degree (o)','Fontsize',12);
ylabel('halfpeakTime (s)','Fontsize',12);
grid on;
 subplot(2,1,2);
 plot(x,halfTime2,'-*');xlabel('Degree (o)','Fontsize',12);
ylabel('halfpeakTime (s)','Fontsize',12); title('datatiffing')
grid on;
%%
figure(16);hold off;grid on;
data2Save=zeros(p,px2+1);
for ii=1:px2
    ii2=centerP3(ii);
    dataPlot=data2(:,ii2);
    dataPlot=dataPlot+minAll;
    if ii==2
        hold on;
    end
        plot(timeCourse,dataPlot)
        plot(timeCourse(Ind(ii2)),dataPlot(Ind(ii2)),'o')

    minAll=minAll+(dataMin(ii2))*1.5;
%     minAll=minAll+0.08;
%     minAll=minAll+0.03;
    data2Save(:,ii+1)=dataPlot;
end
data2Save(:,1)=timeCourse(:);
title('IOS curve')

figure(17);hold off; subplot(2,1,1); plot(x,abs(dataMin),'-*');xlabel('Degree (o)','Fontsize',12);
ylabel('IOS (\DeltaI/I)','Fontsize',12);
grid on;
title('IOS_max')
subplot(2,1,2);

plot(x,abs(maxValueFitting),'-*');xlabel('Degree (o)','Fontsize',12);
ylabel('IOS (\DeltaI/I)','Fontsize',12);
grid on;
title('IOS_datafitting')
peakTime=timeCourse(Ind);
% figure(20);
% % subplot(2,1,1); 
% plot(x,peakTime,'-*');xlabel('Degree (o)','Fontsize',12);
% ylabel('peakTime (s)','Fontsize',12);


%% dealing with half peak


%% get statistics
data3=data;
for ii=1:px
    data3(:,ii)=(data3(:,ii)-mean(data3(1:pre,ii)));
end
[dataMin3,Ind]=min(data3(searZone,:));
figure(18); hold off; plot(x,abs(dataMin3),'-*');
xlabel('Degree (o)','Fontsize',12);
ylabel('absolute IOS (\DeltaI)','Fontsize',12);
grid on;
title('abslute')
result=['angle_',num2str(interval)];
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end
saveas(19,fullfile(resultPath,'halfPeaktime.tif'),'tiff')
saveas(18,fullfile(resultPath,'absolute.tif'),'tiff')
saveas(17,fullfile(resultPath,'IOS.tif'),'tiff')
saveas(16,fullfile(resultPath,'IOSraw.tif'),'tiff')
saveas(15,fullfile(resultPath,'Imgraw.tif'),'tiff')
saveas(4,fullfile(resultPath,'datafitting.tif'),'tiff')
save(fullfile(resultPath,'absolute.txt'),'dataMin3','-ASCII')
dataMinSave=[x.',abs(dataMin(:))];
save(fullfile(resultPath,'IOS_detaI_I.txt'),'dataMinSave','-ASCII')

dataMinSave=[x.',abs(maxValueFitting(:))];
save(fullfile(resultPath,'IOS_detaI_I_fitting.txt'),'maxValueFitting','-ASCII')

save(fullfile(resultPath,'halfTime.txt'),'halfTime','-ASCII')
save(fullfile(resultPath,'halfTime_byfittingExpo.txt'),'halfTime2','-ASCII')
dataMin3Save=[x.',abs(dataMin3(:))];
save(fullfile(resultPath,'IOS_absolute.txt'),'dataMin3Save','-ASCII')
save(fullfile(resultPath,'IOS_raw.txt'),'data','-ASCII')
save(fullfile(resultPath,'IOS_fewRegion.txt'),'data2Save','-ASCII')
function [timeToHalf,varargout]=findHalf_expo(x0,y0,n1,n2,varargin)
baseY=mean(y0(n1-99:n1));
if abs(x0(n1))>0.0001
    error('time scale is not right')
end
x=x0(n1:n2);
x=x(:);
y=y0(n1:n2);
y=y(:);
fun=inline('y-baseY+par(1)-par(1)*exp(-x/par(2))','par','y','x','baseY');
fun2=inline('baseY-par(1)+par(1)*exp(-x/par(2))','par','x','baseY');
par0=[0.01    0.001];
[par,resnorm,residual]=lsqnonlin(fun,par0,[],[],[],y,x,baseY);
% figure(3);hold off; plot(x,fun2(par,x,baseY),x,y,'*')
if ~isempty(varargin)
    y2=fun2(par,x,baseY);
    y3=fun2(par,zeros(n1-1,1),baseY);
    y4=y0(1:n2);y4(1:n1-1)=y3;y4(n1:end)=y2;
    
    figure(4);hold off; plot(x0,y0,'+');grid on;
    hold on; plot(x0(1:n2),y4,'r-','lineWidth',2)
end
timeToHalf=-par(2)*log(0.5);
if nargout>1
    maxValue=baseY-par(1);
    varargout{1}=maxValue;
    
end


function me_continueShow_Callback(hObject, eventdata, handles)
continousROIsShow(handles);


function pushbutton131_Callback(hObject, eventdata, handles)
handles=pushbutton31_Callback(hObject, eventdata, handles);
handles=pushbutton32_Callback(hObject, eventdata, handles);
handles=pushbutton81_Callback(hObject, eventdata, handles);
handles=pushbutton79_Callback(hObject, eventdata, handles);
handles=pushbutton80_Callback(hObject, eventdata, handles);


function listbox3_Callback(hObject, eventdata, handles)


function listbox3_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton132_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.*','select movie please');
cd(filePath)
aviObject=VideoReader(file);
p=aviObject.NumberOfFrames;
dotNO=strfind(file,'.');
dotNO=dotNO(end);
fileType=file(dotNO+1:end);
fileTmp=file(1:dotNO-1);
%% result Path
result=fileTmp;
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end
%% see the first frame
image1=read(aviObject,1);
figure(3); imshow(image1,[]);
[image1Ind,colorMap_1]=rgb2ind(image1,256);
[mm,nn]=size(image1(:,:,1));
x1=1;y1=1;
x2=nn;y2=mm;
% x1=288; y1=1;
% x2=930; y2=835;
% x1=120; y1=189;
% x2=1082;y2=800;
% x1=231; x2=894;y1=1;y2=mm;
h_wait=waitbar(0,'please wait');
for ii=1:p
    waitbar(ii/p,h_wait,[num2str(ii*100/p,'%03.1f'),'completed']);
    image1=read(aviObject,ii);
    imwrite(image1(:,:,1),fullfile(resultPath,[fileTmp,num2str(ii,'%03d'),'.tif']),'tif','compression','none')
%     image1Ind=rgb2ind(image1,colorMap_1);
%     X=image1Ind(y1:y2,x1:x2);
%     MAP=colorMap_1;
%      M(ii)=im2frame(X,MAP);
end
close(h_wait)


function me_3d_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
result=['threeD'];
resultPath=fullfile(filePath,result);
if exist(resultPath)==7 
   
else
    mkdir(resultPath)
end
b=differentTypeRead(fullfile(filePath,file),handles.fileType);
b=double(b);
% b=b(end:-1:1,:);
% b(1:20,1:20)=255;
[m,n]=size(b);
x=1:n;
y=1:m;
[xx,yy]=meshgrid(x,y);
xi=1:.2:n;
yi=1:.2:m;
[xxi,yyi]=meshgrid(xi,yi);
zi=interp2(xx,yy,b,xxi,yyi);
% zi = imresize(zi,0.1,'box');

% figure; surf(zi)
% colorbar

xt=200;
kernel=fspecial('average',[xt xt]);
zi_2=zi;
xt=100
[mm,nn]=size(zi_2);
kernel=fspecial('gaussian',5,3);
% zi_2=imfilter(zi_2,kernel);
%% replicate
zi_2=imfilter(zi_2,kernel,'replicate');
% zi_2=imfilter(zi_2,kernel,'replicate');
% for jj=1:mm
%     for ii=1:nn
%        if zi_2(jj,ii)>240
%             zi_2(jj,ii)=zi_2(jj,ii)-50;
%         end
%     end
% end 

% figure;imshow(uint8(zi_2));
% 
% h.fig=figure; h.surf=surf(zi_2);colorbar
% set(gca,'zlim',[-80,30000])
% 
% h.axis=get(h.surf,'parent');
% 
% 
% set(h.axis,'box','on')
% set(h.axis,'color',[.8,.8,.8])
% 
% 
% h2=h;
% 
% 
h.fig=figure(1);hold off;
set(gca,'unit','pixel')
%% crop
% a=3;
% h.surf=mesh(zi_2(a:end-a,a:end-a));h.h=colorbar;
h.surf=mesh(zi_2);h.h=colorbar;
h.axis=get(h.surf,'parent');
set(gcf,'color',[1 1 1])
set(gca,'zlim',[min(b(:)),max(b(:))])
set(h.axis,'Visible','off')
% colormap('gray')
% grid off
set(h.axis,'color',[1 1 1])
set(h.surf,'facecolor','interp')
set(h.surf,'EdgeColor','none','FaceLighting','phong')
set(h.surf,'edgeAlpha',0)
% daspect([5 5 1])
% axis tight
% view(-50,30)
% camlight left
%%
% set(h.axis,'clim',[108,148])
% set(h.fig,'Color',[.6,.6,.6])
% a=(-0.20:0.05:0.25);
% set(h.h,'YTickLabel',a)
az =-5
el =30
view(az, el);
% set(h.axis,'visible','off')
cd(resultPath);
colorbar off
% saveas(gcf,'threeD','tiff')
% saveas(gcf,'threeD','bmp')
% print('-dtiff',['three_dimension_',file])
% cd ..
% export_fig test.png -transparent


function cb_ms_dI_Callback(hObject, eventdata, handles)


function cb_pnIOS_filter_Callback(hObject, eventdata, handles)


function cb_IOS_filter_Callback(hObject, eventdata, handles)


function Untitled_20_Callback(hObject, eventdata, handles)



function pushbutton133_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
shift=handles.shift;
[resultFile,resultPath]=img_register(file,filePath,handles.fileName,shift,'bilinear',handles);
function [resultFile,resultPath]=img_register(file,filePath,fileName,shift,method,handles)
subShift=shift;
subShift=subShift-floor(subShift);

shift=floor(shift);
shift=-shift;

cd(filePath)
result='registeredImages';
resultPath=fullfile(filePath,result);
cd(filePath)
if exist(resultPath)==7
else
    mkdir(result)
end
dotNO=strfind(file,'.');
fileType=file(dotNO(end)+1:end);
% fileName=dir(['*.',fileType]);for iss=length(fileName):-1:1;if strcmp(fileName(iss).name(1:2),'._'); fileName(iss)=[];end;end
p=length(fileName);
img2=imread(file);
img=img2(:,:,1);
imgHandle=imfinfo(file);
BitDepth=imgHandle.BitDepth;

minShift=abs(min(shift));
maxShift=max(shift);

yPlus=max(maxShift(1),minShift(1));
xPlus=max(maxShift(2),minShift(2));

[m,n]=size(img);
h_wait=waitbar(0,'please wait');
fileType=handles.fileType;
imgDepth=handles.imgDepth;
% filePath
for ii=1:p
    waitbar(ii/p,h_wait,[num2str(ii/p*100,'%04.1f'),'%completed']);
    file=fileName(ii).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);

    yStart=shift(ii,1)+yPlus+1;
    xStart=shift(ii,2)+xPlus+1;
    yEnd=shift(ii,1)+yPlus+m;
    xEnd=shift(ii,2)+xPlus+n;
    imgTmp=zeros(m+2*yPlus,n+2*xPlus);
%    imgTmp=zeros(m+2*yPlus,n+2*xPlus,['uint',num2str(BitDepth)]);

    imgTmp(yStart:yEnd,xStart:xEnd)=img;
    imgTmp2=img_ROI_registerSubpixel(imgTmp,subShift(ii,:),method);
    differentTypeWrite(imgTmp2,fullfile(resultPath,file),fileType,imgDepth);
%     imgTmp3=eval(['uint',num2str(BitDepth),'(imgTmp2)']);
%     cd(resultPath)
%     dotNO=strfind(file,'.');
%     saveFile=file(1:dotNO(end));
%     imwrite(imgTmp3,[saveFile,'tif'],'tiff','compression','none');
    
end
close(h_wait)
cd(resultPath)
shiftSave=zeros(p,3);
shiftSave(:,1)=(1:p).';
shiftSave(:,2:end)=shift;
save('shift.txt','shiftSave','-ASCII')
cd ..


function handles=pushbutton134_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
ref=fft2(img);
result='fftRegis';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end

h_wait=waitbar(0,'please wait');
a=eval(get(handles.et_fftRegsi,'string'));
shift=zeros(p,2);
for ii=1:p
    waitbar(ii/p,h_wait,[num2str(ii/p*100,'%04.1f'),'%completed']);
    file=fileName(ii).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    [output, img2] = dftregistrationIOS(ref,fft2(double(img)),a);
    shift(ii,:)=output(3:4);
    differentTypeWrite(abs(ifft2(img2)),fullfile(resultPath,file),fileType,handles.imgDepth);
end

close(h_wait)
figure;plot(1:p,shift(:,1),'r',1:p,shift(:,2),'b');
legend('shift on y','shift on x')
% xlswrite(fullfile(resultPath,'shift_y_x.xls'),shift);
% save(fullfile(resultPath,'shift_y_x.txt'),'shift','-ASCII');


function [output Greg] = dftregistrationIOS(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk 
% and James R. Fienup. 
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued 
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft))); 
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero); 
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    output=[error,diffphase];
        
% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc); 
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    md2 = fix(m/2); 
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end

    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
% Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
  
    % Compute crosscorrelation and locate the peak 
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak 
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2 
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;

    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;     
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftupsIOS(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid 
        [max1,loc1] = max(CC);   
        [max2,loc2] = max(max1); 
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftupsIOS(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftupsIOS(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);  
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    

    % If upsampling = 2, no additional pixel shift refinement
    else    
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end  

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end


function out=dftupsIOS(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;


function et_fftRegsi_Callback(hObject, eventdata, handles)


function et_fftRegsi_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cb_cdCurrent_Callback(hObject, eventdata, handles)


function Untitled_21_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.*','get background img');
dotNO=strfind(file,'.');fileType0=file(dotNO(end)+1:end);
img0=differentTypeRead(fullfile(filePath,file),fileType0);
handles.minusImg0=img0;
guidata(hObject,handles)

function Untitled_22_Callback(hObject, eventdata, handles)

img0=handles.minusImg0;
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
result='minusImg';
resultFolder=fullfile(handles.filePath,result);

if exist(resultFolder)==7
else
    cd(filePath)
    mkdir(result)
end
cd(resultFolder)
[fileTmp,resultFolder]=uiputfile(file);

h_wait=waitbar(0,'wait');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
%     cd(filePath)
    img=differentTypeRead(fullfile(filePath,file),fileType);
    vv=1;

    differentTypeWrite(img-img0,file,fileType,handles.imgDepth);
end
close(h_wait)

guidata(hObject,handles)


function et_IOSbigger_Callback(hObject, eventdata, handles)


function et_IOSbigger_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cb_biggerBack_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.*','select movie please');
cd(filePath)
aviObject=VideoReader(file);
p=aviObject.NumberOfFrames;
dotNO=strfind(file,'.');
dotNO=dotNO(end);
fileType=file(dotNO+1:end);
fileTmp=file(1:dotNO-1);
%% see the first frame
image1=read(aviObject,1);
figure(3); imshow(image1,[]);
[image1Ind,colorMap_1]=rgb2ind(image1,256);
[mm,nn]=size(image1(:,:,1));
x1=1;y1=1;
x2=nn;y2=mm;

result='imgs';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end


% x1=288; y1=1;
% x2=930; y2=835;
% x1=120; y1=189;
% x2=1082;y2=800;
% x1=231; x2=894;y1=1;y2=mm;
h_wait=waitbar(0,'please wait');
for ii=1:p
    waitbar(ii/p,h_wait,[num2str(ii*100/p,'%03.1f'),'completed']);
    image1=read(aviObject,ii);
%     save()
%     image1Ind=rgb2ind(image1,colorMap_1);
%     X=image1Ind(y1:y2,x1:x2);
%     MAP=colorMap_1;
%      M(ii)=im2frame(X,MAP);
     imwrite(image1,fullfile(resultPath,['img',num2str(ii,'%03d'),'.tif']),'tif');
end
close(h_wait)


% --------------------------------------------------------------------
function Untitled_24_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Untitled_25_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
% cd(filePath)
img=differentTypeRead(fullfile(filePath,file),fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
rawData_x=zeros(p,sum(ROI_N(:))+1);

timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData_x(:,1)=timeCourse;
rawData_y=rawData_x;
rawData_ratio=rawData_x;
mean_x=rawData_x;
mean_y=rawData_x;
h_wait=waitbar(0,'wait');

x1=1:n;
y1=1:m;
[xx,yy]=meshgrid(x1,y1);


for kk=1:p
    %h_wait=waitbar(0,'wait');
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
 if get(handles.pp_spfilterMethod,'value')>1
    img=differentTypeReadFilter_handles(img,handles);
 end     
    vv=1;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
            vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bw=roipoly(img(:,:,1),cx,cy);
            ROItmp=std(img(bw)/sum(img(bw)).*xx(bw));
            rawData_x(kk,vv)=ROItmp;
            mean_x(kk,vv)=sum(img(bw)/sum(img(bw)).*xx(bw))/sum(bw(:));
            ROItmp=std(img(bw)/sum(img(bw)).*yy(bw));
            rawData_y(kk,vv)=ROItmp;   
            mean_y(kk,vv)=sum(img(bw)/sum(img(bw)).*yy(bw))/sum(bw(:));
%             rawData_ratio(kk,vv)=mean_y(kk,vv)/rawData_x(kk,vv);
        end
    end
end
close(h_wait)
rawData_x(:,2:end)= ...
    IOS_time_gui_filter(rawData_x(:,2:end),handles,'vertical');
rawData_y(:,2:end)= ...
    IOS_time_gui_filter(rawData_y(:,2:end),handles,'vertical');
rawData_ratio(:,2:end)=rawData_y(:,2:end)./rawData_x(:,2:end);
figure(18);close(18);figure(18);
subplot(1,3,1);plot(rawData_ratio(:,1),rawData_x(:,2:end));title('SD of x')
subplot(1,3,2);plot(rawData_ratio(:,1),rawData_y(:,2:end));title('SD of y')
subplot(1,3,3);plot(rawData_ratio(:,1),rawData_ratio(:,2:end));title('SD of y/SD of x')
% rawData=log(rawData);
% handles.rawData_meanOfROI=rawData;
guidata(hObject,handles)


function cb_circle_Callback(hObject, eventdata, handles)


function cb_gr11_Callback(hObject, eventdata, handles)


function cb_circle2_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function pushbutton55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function handles=pb_oneButton_Callback(hObject, eventdata, handles)
set(handles.pb_oneButton,'backgroundcolor',handles.bgColor);
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
ref=fft2(img);
result='fftRegis';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end

h_wait=waitbar(0,'Wait','Name','Registration in progress');
a=eval(get(handles.et_fftRegsi,'string'));
shift=zeros(p,2);
for ii=1:p
    waitbar(ii/p,h_wait,[num2str(ii/p*100,'%04.1f'),'%completed']);
    file=fileName(ii).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    img=differentTypeReadFilter_handles(img,handles);
    [output, img2] = dftregistrationIOS(ref,fft2(double(img)),a);
    shift(ii,:)=output(3:4);
    differentTypeWrite(abs(ifft2(img2)),fullfile(resultPath,file),fileType,handles.imgDepth);
end

close(h_wait)
% figure;plot(1:p,shift(:,1),'r',1:p,shift(:,2),'b');
% legend('shift on y','shift on x')
% xlswrite(fullfile(resultPath,'shift_y_x.xls'),shift);
% save(fullfile(resultPath,'shift_y_x.txt'),'shift','-ASCII');

handles.filePath_fftRegis=resultPath;
handles.filePath = resultPath;

handles=pb_rawIOS_mat_Callback(hObject, eventdata, handles);
handles.filePath=handles.mat_resultPath;
handles.fileType='mat';
handles.fileName=dir(fullfile(handles.filePath,'*.mat'));
handles.file=handles.fileName(1).name;
handles.p=length(handles.fileName);

handles.filePathMat=handles.mat_resultPath;
handles.fileNameMat=dir(fullfile(handles.filePath,'*.mat'));
handles.fileMat=handles.fileName(1).name;
handles.pMat=length(handles.fileNameMat);

p=handles.p;
listName=cell(p,1);
for ii=1:p
    listName{ii}=handles.fileName(ii).name;
    
end
set(handles.dt_IntestedImg,'string',['1:',num2str(p)]);
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
set(handles.pb_oneButton,'BackgroundColor',[1,0.49,0.28]);
handles=imageShow(handles);
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles)

function handles=handlesUpdata(handles,eventdata)
handles.file=handles.fileName(1).name;
handles.p=length(handles.fileName);
p=handles.p;
listName=cell(p,1);
for ii=1:p
    listName{ii}=handles.fileName(ii).name;
    
end
set(handles.dt_IntestedImg,'string',['1:',num2str(p)]);
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
handles=imageShow(handles);

handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);

function handles=pb_Combined_images_Callback(hObject, eventdata, handles)
file=handles.fileMat;
filePath=handles.filePathMat;
fileName=handles.fileNameMat;
p=length(fileName);
fileType='mat';%handles.fileType;

% set(handles.dt_IntestedImg,'string',['1:',num2str(p)])
% set(handles.lb_curFileNames,'string',fileName)
% set(handles.lb_curFileNames,'value',1)
% handles=imageShow(handles);


img=differentTypeRead(fullfile(filePath,file),fileType);
% ref=fft2(img);
result='Combined_images';
resultPath=fullfile(handles.filePath_output,result);
if exist(resultPath)~=7
    mkdir(resultPath)
end

if exist(handles.filePath_fftRegis)
    filePath1=[handles.filePath_fftRegis,filesep];
    fileName1=dir([filePath1,'*.tif']);
else
    filePath1=[handles.filePath_output,filesep];
    fileName1=dir([filePath1,'*.tif']);
    if exist(handles.oriRange)
        fileName1 = fileName1(handles.oriRange);
    end
end



% result1='grayBackgroundImages';
% resultPath1=fullfile(filePath1,result1);
% mkdir(resultPath1);

% aa=differentTypeRead(fullfile(filePath1,file),fileType);
aa=imread([filePath1,fileName1(1).name]);%fileName(1).name = a00000.BMP;
aa=aa(:,:,1);
[m,n]=size(aa);
img2=uint8(zeros(m,n));


p2=size(fileName1,1);
loop=floor(p2/p);

%
h_wait=waitbar(0,'please wait');
%


D=uint8(zeros(m,n,3));
B=uint8(zeros(m,n,3));
ImgTT=uint8(zeros(m,n,3));
D(:,:,1)=aa;
D(:,:,2)=aa;
D(:,:,3)=aa;

    cd(resultPath);
for tt=1:p
    waitbar(tt/p,h_wait,[num2str(tt/p*100,'%04.1f'),'%completed']);
    imgTemp=double(zeros(m,n));
    for ii=((tt-1)*loop+1):(tt*loop)
%     for ii= tt:p:p2
        file1=fileName1(ii).name;
%         fileName1(ii).name
        img1=differentTypeRead(fullfile(filePath1,file1),'tif');
%         img1=differentTypeRead(fullfile(filePath1,file1),'mat');
        %     img=differentTypeReadFilter_handles(img,handles);
        imgTemp=double(imgTemp+img1);
    end
    img2=im2uint8(mat2gray(imgTemp));
    
    file=fileName(tt).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    newIntensity=inline(get(handles.et_inline,'string'));
    contents = cellstr(get(handles.pp_IOScolormap,'String'));
    colorSelectedTmp=contents(get(handles.pp_IOScolormap,'value'));
    colorSelected=eval(colorSelectedTmp{1});
    img3=newIntensity(img);
    imgColor=ind2rgb(uint8(img3),colorSelected);
    
    
    
    B=im2uint8(imgColor); % Background
 
    D(:,:,1)=img2; % Difference
    D(:,:,2)=img2;
    D(:,:,3)=img2;
    
    for ss=1:m
        for zz=1:n
            if (B(ss,zz,1)==0 && B(ss,zz,2)<255)||(B(ss,zz,3)==0 && B(ss,zz,2)<255)  %%%blue||red
                ImgTT(ss,zz,:)=B(ss,zz,:);
            else
                ImgTT(ss,zz,:)=D(ss,zz,:);
                
            end
        end
    end
    

    indexNO=strfind(file1,'.tif');
    
    scaleName=get(handles.et_inline,'string');
    dotNO1=strfind(scaleName,'/');
    scale=scaleName(dotNO1+1:end);
%     imwrite(ImgTT,['combined_',file1(1:indexNO(end)-1),'.tif']);
 imwrite(ImgTT,['combined_scale',scale,'_',file1(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif']);
end

handles.filePath=resultPath;
handles.fileType='tif';
cd(resultPath);
fileName=dir(['*.','tif']);
handles.fileName=fileName;
handles.file=['combined_scale',scale,'_',file1(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif'];
file=['combined_scale',scale,'_',file1(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif'];
handles.p=length(handles.fileName);
p=handles.p;

listName=cell(p,1); 

for ii=1:p
    listName{ii}=handles.fileName(ii).name;
end

set(handles.dt_IntestedImg,'string',['1:',num2str(p)]);
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
handles=imageShow(handles);
close(h_wait);
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles);
 
function pb_Combined_NoBlue_Callback(hObject, eventdata, handles)
file=handles.fileMat;
filePath=handles.filePathMat;
fileName=handles.fileNameMat;
p=length(fileName);
fileType='mat';%handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
% ref=fft2(img);
result='Combined_images_NoBlue';
resultPath=fullfile(handles.filePath_output,result);
if exist(resultPath)~=7
    mkdir(resultPath)
end


if exist(handles.filePath_fftRegis)
    filePath1=[handles.filePath_fftRegis,filesep];
    fileName1=dir([filePath1,'*.tif']);
else
    filePath1=[handles.filePath_output,filesep];
    fileName1=dir([filePath1,'*.tif']);
    if exist(handles.oriRange)
        fileName1 = fileName1(handles.oriRange);
    end
end
% result1='grayBackgroundImages';
% resultPath1=fullfile(filePath1,result1);
% mkdir(resultPath1);
aa=imread([filePath1,fileName1(1).name]);%fileName(1).name = a00000.BMP;
aa=aa(:,:,1);
[m,n]=size(aa);
img2=uint8(zeros(m,n));
imgTemp=double(zeros(m,n));

p2=size(fileName1,1);
loop=floor(p2/p);

%
h_wait=waitbar(0,'please wait');
%


D=uint8(zeros(m,n,3));
B=uint8(zeros(m,n,3));
ImgTT=uint8(zeros(m,n,3));
D(:,:,1)=aa;
D(:,:,2)=aa;
D(:,:,3)=aa;

cd(resultPath);
for tt=1:p
    waitbar(tt/p,h_wait,[num2str(tt/p*100,'%04.1f'),'%completed']);
    for ii=((tt-1)*loop+1):(tt*loop)
        
        file1=fileName1(ii).name;
        img1=differentTypeRead(fullfile(filePath1,file1),'tif');
        %     img=differentTypeReadFilter_handles(img,handles);
        imgTemp=double(imgTemp+img1);
    end
    img2=im2uint8(mat2gray(img1));
    
    file=fileName(tt).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    newIntensity=inline(get(handles.et_inline,'string'));
    contents = cellstr(get(handles.pp_IOScolormap,'String'));
    colorSelectedTmp=contents(get(handles.pp_IOScolormap,'value'));
    colorSelected=eval(colorSelectedTmp{1});
    img3=newIntensity(img);
    imgColor=ind2rgb(uint8(img3),colorSelected);
    
    
    
    B=im2uint8(imgColor);
    
    D(:,:,1)=img2;
    D(:,:,2)=img2;
    D(:,:,3)=img2;
    for ss=1:m
        for zz=1:n
            if (B(ss,zz,3)==0 && B(ss,zz,2)<255)  %%%blue||red
                ImgTT(ss,zz,:)=B(ss,zz,:);
            else
                ImgTT(ss,zz,:)=D(ss,zz,:);
                
            end
        end
    end
    
    indexNO=strfind(file1,'.');
    scaleName=get(handles.et_inline,'string');
    dotNO1=strfind(scaleName,'/');
    scale=scaleName(dotNO1+1:end);
    imwrite(ImgTT,['combined-noblue_scale',scale,'_',file1(1:indexNO(end)-1),'.tif']);
end

handles.filePath=resultPath;
handles.fileType='tif';
cd(resultPath);
fileName=dir(['*.','tif']);
handles.fileName=fileName;
handles.file=['combined-noblue_scale',scale,'_',file1(1:indexNO(end)-1),'.tif'];
file=['combined-noblue_scale',scale,'_',file1(1:indexNO(end)-1),'.tif'];
handles.p=length(handles.fileName);
p=handles.p;

listName=cell(p,1);

for ii=1:p
    listName{ii}=handles.fileName(ii).name;
end
set(handles.dt_IntestedImg,'string',['1:',num2str(p)]);
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
handles=imageShow(handles);
close(h_wait);
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles)
 
% hObject    handle to pb_Combined_NoBlue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


 
% hObject    handle to pb_Combined_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function pushbutton136_Callback(hObject, eventdata, handles)
[file,filePath]=uiputfile('*.mat');
save(fullfile(filePath,file),'handles')

function pushbutton137_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.*');
h=load(fullfile(filePath,file));
guidata(hObject,handles)



function edit109_Callback(hObject, eventdata, handles)
% hObject    handle to edit109 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit109 as text
%        str2double(get(hObject,'String')) returns contents of edit109 as a double


% --- Executes during object creation, after setting all properties.
function edit109_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit109 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton138.
function pushbutton138_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton138 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit110_Callback(hObject, eventdata, handles)
% hObject    handle to edit110 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit110 as text
%        str2double(get(hObject,'String')) returns contents of edit110 as a double


% --- Executes during object creation, after setting all properties.
function edit110_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit110 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton139.
function pushbutton139_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton139 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_26_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_Combined_images.


% --- Executes on button press in pb_Combined_NoBlue.


% --- Executes on button press in pushbutton145.
function pushbutton145_Callback(hObject, eventdata, handles)
file=handles.file;
filePath_O=handles.filePath;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);

fileName1=dir('*.tif');


aa=imread(fileName1(1).name);%fileName(1).name = a00000.BMP;
aa=aa(:,:,1);
[m,n]=size(aa);
img2=uint8(zeros(m,n));
imgTemp=double(zeros(m,n));

p2=size(fileName1,1);

result=['ratio images'];
mkdir(result)
resultPath=fullfile(filePath,result);
cd(resultPath);
cd(filePath);

C=double(zeros(m,n));
Y=double(zeros(m,n));
Img=double(zeros(m,n));


for kk=1:p/2
    cd(filePath);
    C(:,:)=imread(fileName(kk).name);
    Y(:,:)=imread(fileName(p/2+kk).name);
    Img(:,:)=Y./C;
    cd(resultPath);
%     imwrite(Img,fileName(kk).name);
    saveName=[fileName(kk).name(1:end-3),'mat'];
    save(saveName,'Img');
end
guidata(hObject,handles)
% handles.filePath=resultPath;
% 
% file=handles.file;
% filePath=handles.filePath;
% fileName=handles.fileName;
% p=length(fileName);
% fileType=handles.fileType;
% img=differentTypeRead(fullfile(filePath,file),fileType);
% ref=fft2(img);
% result2='fftRegis';
% resultPath2=fullfile(filePath,result2);
% if exist(resultPath2)==7
% else
%     mkdir(resultPath2)
% end
% 
% h_wait=waitbar(0,'please wait');
% a=eval(get(handles.et_fftRegsi,'string'));
% shift=zeros(p,2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imagN=str2num(get(handles.et_imageNo,'string'));
% for ii=(p2/2-imagN+1):p2/2
%     waitbar(ii/p,h_wait,[num2str(ii/p*100,'%04.1f'),'%completed']);
%     cd(resultPath);
%     file=fileName(ii).name;
%     img=differentTypeRead(fullfile(filePath,file),fileType);
%     img=differentTypeReadFilter_handles(img,handles);
%     [output, img2] = dftregistrationIOS(ref,fft2(double(img)),a);
%     shift(ii,:)=output(3:4);
%     differentTypeWrite(abs(ifft2(img2)),fullfile(resultPath2,file),fileType,handles.imgDepth);
% end
% 
% close(h_wait)
% %%very important:change everything in the handle before you use the pb_rawIOS_mat
% 
% handles.filePath=resultPath2;
% handles.fileName=dir(fullfile(handles.filePath,'*.tif'));
% handles.file=handles.fileName(1).name;
% fileName=handles.fileName;
% p=length(fileName);
% handles.p=p;
% 
% handles=pb_rawIOS_mat_Callback(hObject, eventdata, handles);
% 
% % guidata(hObject,handles)
% 
% for kk=(p2/2-imagN+1):p2/2    
%     cd(filePath_O);
%     handles.fileName=dir(fullfile(filePath_O,'*.tif'));
%     fileName=handles.fileName;
%     Img(:,:)=imread(fileName(kk).name);
% %     Img=im2uint8(mat2gray(double(Img)));
%     cd(resultPath2);
%     imwrite(Img,fileName(kk).name);
% end
% 
% handles.filePath=handles.mat_resultPath;
% handles.fileType='mat';
% handles.fileName=dir(fullfile(handles.filePath,'*.mat'));
% handles.file=handles.fileName(1).name;
% handles.p=length(handles.fileName);
% p=handles.p;
% listName=cell(p,1);
% for ii=1:p
%     listName{ii}=handles.fileName(ii).name;
%     
% end
% set(handles.dt_IntestedImg,'string',['1:',num2str(p)])
% set(handles.lb_curFileNames,'string',listName)
% set(handles.lb_curFileNames,'value',1)
% set(handles.st_curFilePath,'string',handles.filePath);
% handles=imageShow(handles);
% handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
% cd(handles.filePath)
% guidata(hObject,handles)
% hObject    handle to pushbutton145 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function et_imageNo_Callback(hObject, eventdata, handles)
% hObject    handle to et_imageNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_imageNo as text
%        str2double(get(hObject,'String')) returns contents of et_imageNo as a double


% --- Executes during object creation, after setting all properties.
function et_imageNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_imageNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton146.
function pushbutton146_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton146 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file=handles.file;
filePath_O=handles.filePath;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);

fileName1=dir('*.tif');


aa=imread(fileName1(1).name);%fileName(1).name = a00000.BMP;
aa=aa(:,:,1);
[m,n]=size(aa);
img2=uint8(zeros(m,n));
imgTemp=double(zeros(m,n));

p2=size(fileName1,1);

result1=['ratio images_bgs'];
result2=['cfp_bgs'];
result3=['yfp_bgs'];
mkdir(result1)
resultPath1=fullfile(filePath,result1);
mkdir(result2)
resultPath2=fullfile(filePath,result2);
mkdir(result3)
resultPath3=fullfile(filePath,result3);
cd(resultPath1);
cd(filePath);

C=double(zeros(m,n));
C1=uint16(zeros(m,n));
Y=double(zeros(m,n));
Y1=uint16(zeros(m,n));
Img=double(zeros(m,n));


for kk=1:p/2
    cd(filePath);
    C(:,:)=imread(fileName(kk).name);
    C1(:,:)= C(:,:)-min(min(C))+1;
    Y(:,:)=imread(fileName(p/2+kk).name);
    Y1(:,:)= Y(:,:)-min(min(Y))+1;
    Img(:,:)=Y1./C1;
    cd(resultPath1);
%     imwrite(Img,fileName(kk).name);
    saveName=[fileName(kk).name(1:end-3),'mat'];
    save(saveName,'Img');
    cd(resultPath2);
%     saveName=[fileName(kk).name(1:end-3),'mat'];
%     save(saveName,'C');
    imwrite(C1,fileName(kk).name);
    cd(resultPath3);
    imwrite(Y1,fileName(kk).name);
%     saveName=[fileName(kk).name(1:end-3),'mat'];
%     save(saveName,'Y');
end
guidata(hObject,handles)




% --- Executes on button press in pushbutton147.

% hObject    handle to pushbutton147 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function pushbutton147_Callback(hObject, eventdata, handles)
[file,filePath]=uigetfile('*.mat');
load(fullfile(filePath,file))
handles.rawData_positive_meanOfROI=rawData_positive_meanOfROI;
handles.rawData_pn_meanOfROI=rawData_pn_meanOfROI;
handles.rawData_negative_meanOfROI=rawData_negative_meanOfROI;
handles.pointPst=pointPst;
handles.positiveMap=positiveMap;
handles.negativeMap=negativeMap;
handles.colorTmp=colorTmp;
guidata(hObject,handles)

% --- Executes on button press in pushbutton148.
function handles=pushbutton148_Callback(hObject, eventdata, handles)
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
pnYellow=and(positiveMap,handles.negativeMap);
positiveMap2=xor(pnYellow,positiveMap);
negativeMap2=xor(pnYellow,negativeMap);
%% positive
rawData=handles.rawData_positive_meanOfROI;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;

legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
%% draw all together
if size(rawData,2)==2
    pre=eval(get(handles.et_pre2,'string'));
    meanP=mean(rawData(1:pre,2));

    pPlot=(rawData(:,2)-meanP)/meanP;
    
    nPlot=handles.rawData_negative_meanOfROI(:,2);
    meanN=mean(nPlot(1:pre));
    nPlot=(nPlot-meanN)/meanN;
    
    
    cx=pointPst{1}{1}(:,1);
    cy=pointPst{1}{1}(:,2);
    bwTmp=roipoly(double(positiveMap),cx,cy);
    bwTmpNumber=and(bwTmp,positiveMap2);
    pNumber=sum(bwTmpNumber(:));
    bwTmpNumber=and(bwTmp,negativeMap2);
    nNumber=sum(bwTmpNumber(:));
    pnPlot=(rawData(:,2)*pNumber+handles.rawData_negative_meanOfROI(:,2)*nNumber)/(nNumber+pNumber);
    meanPN=mean(pnPlot(1:pre));
    pnPlot=(pnPlot-meanPN)/meanPN;
    pPlot= IOS_time_gui_filter(pPlot,handles,'vertical');
    nPlot= IOS_time_gui_filter(nPlot,handles,'vertical');
    figure(108);
    plot(timeCourse,pPlot,'r',timeCourse,nPlot,'g','lineWidth',3)
    xlabel(timeCourseUnit)
    handles.rawData_pnPlot=[timeCourse(:),pPlot(:),nPlot(:)];
    
%     plot(timeCourse,pPlot,timeCourse,nPlot,timeCourse,pnPlot)
%     bwTmpNumberSum=sum(bwTmpNumber(:));
end

%% 
pPlot=handles.rawData_positive_meanOfROI(:,2:end);
nPlot=handles.rawData_negative_meanOfROI(:,2:end);
pre=eval(get(handles.et_pre2,'string'));
for ii=1:size(pPlot,2)
    meanP=mean(pPlot(1:pre,ii));
    if abs(meanP)>0.00001
        pPlot(:,ii)=(pPlot(:,ii)-meanP)/meanP;
        
    end
    meanN=mean(nPlot(1:pre,ii));
    if abs(meanN)>0.00001
        nPlot(:,ii)=(nPlot(:,ii)-meanN)/meanN;
        
    end    
end
handles.rawData_pnPlot=[timeCourse(:),pPlot,nPlot];
figure(108);close(108);figure(108)
vv=1;
minN=min(nPlot(:));
maxP=max(pPlot(:));
maxP=max([maxP,abs(minN)]);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        subplot(1,size(pPlot,2),vv)
        vv=vv+1;
%         if vv==2
%             hold on;
%         end
        plot(timeCourse,pPlot(:,vv-1),'r-')
        hold on;
        plot(timeCourse,nPlot(:,vv-1),'b-')
        ylim([-maxP,maxP])
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,positiveMap2);
        legend(['group',num2str(ii),'ROI',num2str(jj)]);
    end
end
% legend(legendName)
%%
vv=1;
figure(105);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,positiveMap2);
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum(bwTmpNumber(:)))];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on; title('positive')
handles.rawData_positive_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;
%% negative
rawDataPositive=rawData;
rawData=handles.rawData_negative_meanOfROI;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(106);
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
pre2=eval(get(handles.et_pre2,'string'));
rawDataAbs=zeros(size(rawData));
rawDataAbs(:,1)=timeCourse;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,negativeMap2);
            
            bwTmp2=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber2=and(bwTmp2,positiveMap2);      
            sum1=sum(bwTmpNumber(:));
            sum2=sum(bwTmpNumber2(:));
            rawDataAbs(:,vv)=-sum1*(rawData(:,vv)-mean(rawData(1:pre2,vv))) ...
                +sum2*(rawDataPositive(:,vv)-mean(rawDataPositive(1:pre2,vv)));
            rawDataAbs(:,vv)=rawDataAbs(:,vv)/(sum1+sum2);
            rawDataAbs(:,vv)=rawDataAbs(:,vv)+(mean(rawData(1:pre2,vv))+ ...
                mean(rawDataPositive(1:pre2,vv)))/2;
            cc2=corrcoef(rawDataPositive(:,vv),rawData(:,vv));
            cc=cc2(1,2);
            
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum(bwTmpNumber(:))),'cc',num2str(cc)];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;title('negative')
handles.rawData_negative_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;

%% absolute
rawData=rawDataAbs;
figure(110);
vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                 cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,negativeMap2);
            
            bwTmp2=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber2=and(bwTmp2,positiveMap2);      
            sum1=sum(bwTmpNumber(:));
            sum2=sum(bwTmpNumber2(:));

        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum1+sum2)];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on;title('absolute')
handles.rawData_abs_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;
%% positive and negative
rawData=handles.rawData_pn_meanOfROI;
p=size(rawData,1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData(:,1)=timeCourse;
pointPst=handles.pointPst;
figure(107);
legendName=cell(size(rawData,2)-1,1);
colorTmp=handles.colorTmp;
ROIii={'-','--',':','-.','*-','*--','*:','*-.'};
vv=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        vv=vv+1;
        if vv>2
            hold on;
        end
        plot(rawData(:,1),rawData(:,vv),ROIii{jj},'color',colorTmp(ii,:))
       %plot(rawData(:,1),rawData(:,vv),'*-','color',colorTmp(ii,:))
                        cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(double(positiveMap),cx,cy);
            bwTmpNumber=and(bwTmp,pnYellow);
        legendName{vv-1}=['group',num2str(ii),'ROI',num2str(jj),'pixels', ...
            num2str(sum(bwTmpNumber(:)))];
    end
end
legend(legendName)
% xLimit=get(gca,'xlim');
% yLimit=get(gca,'ylim');
% pre2=eval(get(handles.et_pre2,'string'));
% line([timeCourse(pre2),timeCourse(pre2)],yLimit,'color','k')
% set(gca,'xlim',xLimit)
% set(gca,'ylim',yLimit)
grid on; title('positive and negative mean')
handles.rawData_pn_meanOfROI=rawData;
xlabel(timeCourseUnit)
hold off;
%% map
% figure(107);
% imshow(handles.positiveMap)
%
% hold on;
% title('positive map')
% colorTmp=handles.colorTmp;
% pointPst=handles.pointPst;
% for ii=1:length(pointPst)
%     for jj=1:length(pointPst{ii})
%         cx=pointPst{ii}{jj}(:,1);
%         cy=pointPst{ii}{jj}(:,2);
%
%         hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
% %             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
% %             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
%     end
% end
% hold off;
%
% figure(108);
% imshow(handles.negativeMap)
% hold on;
% title('negative map')
% colorTmp=handles.colorTmp;
% pointPst=handles.pointPst;
% for ii=1:length(pointPst)
%     for jj=1:length(pointPst{ii})
%         cx=pointPst{ii}{jj}(:,1);
%         cy=pointPst{ii}{jj}(:,2);
%
%         hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
% %             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
% %             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
%     end
% end
% hold off;
figure(109);
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
[mm,nn,nn_trash]=size(negativeMap);
pnMap=zeros(mm,nn,3,'uint8');
pnMapTmp=zeros(mm,nn,'uint8');
pnMapTmp(positiveMap)=255;
pnMap(:,:,1)=pnMapTmp;
pnMapTmp=zeros(mm,nn,'uint8');
pnMapTmp(negativeMap)=255;
pnMap(:,:,2)=pnMapTmp;
if isfield(handles,'file') && isfield(handles,'filePath')
    cd(handles.filePath)
    handles.newIntensity=inline(get(handles.et_inline,'string'));
    fileType=handles.fileType;
    file=handles.file;
    img=differentTypeRead(file,fileType);
    newImg=handles.newIntensity(img);
    newImg=im2uint8(mat2gray(newImg));
%     newImg=uint8(newImg);
    newImg(positiveMap)=0;
    newImg(negativeMap)=0;
    pnMap(:,:,1)=pnMap(:,:,1)+newImg;
    pnMap(:,:,2)=pnMap(:,:,2)+newImg;
    pnMap(:,:,3)=pnMap(:,:,3)+newImg;
end
imshow(pnMap);
handles.pnMap=pnMap;
title('red: positive; green: negative; yellow:both')
hold on;
colorTmp=handles.colorTmp;
pointPst=handles.pointPst;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        hh=plot([cx;cx(1)],[cy;cy(1)],'color',colorTmp(ii,:));
%             handles.pointCurrent_h=[handles.pointCurrent_h;hh];
%             handles.pointCurrent_group=[handles.pointCurrent_group;ii];
    end
end
hold off;
guidata(hObject,handles)

% hObject    handle to pushbutton148 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton149.

% hObject    handle to pushbutton149 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function handles=pushbutton149_Callback(hObject, eventdata, handles)
cd(handles.filePath)
result=get(handles.et_positiveNegativeFolder,'string');
    if get(handles.cb_pnIOS_filter,'value')
        result=[result,handles.filterMethod2.saveName];
    end
resultFolder=fullfile(handles.filePath,result);
if exist(resultFolder)==7
else
    mkdir(result)
end
cd(resultFolder)
% [file,filePath]=uiputfile('meanOfROIStatistics.mat');
% cd(filePath)
file='meanOfROIStatistics.mat';
file=file(1:end-4);
rawData_positive_meanOfROI=handles.rawData_positive_meanOfROI;
rawData_negative_meanOfROI=handles.rawData_negative_meanOfROI;
rawData_abs_meanOfROI=handles.rawData_abs_meanOfROI;
rawData_pn_meanOfROI=handles.rawData_pn_meanOfROI;
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
pointPst=handles.pointPst;
colorTmp=handles.colorTmp;
pnMap=handles.pnMap;
rawData_pnPlot=handles.rawData_pnPlot;
    positiveRawAll=handles.positiveRawAll;
    negativeRawAll= handles.negativeRawAll;
save([file,'_positive.txt'],'rawData_positive_meanOfROI','-ASCII');
save([file,'_negative.txt'],'rawData_negative_meanOfROI','-ASCII');
save([file,'_abs.txt'],'rawData_abs_meanOfROI','-ASCII');
save([file,'pn_debase.txt'],'rawData_pnPlot','-ASCII')
save([file,'_positiveNegative.txt'],'rawData_pn_meanOfROI','-ASCII');
save([file,'.mat'],'rawData_positive_meanOfROI','rawData_negative_meanOfROI', ...
    'positiveMap','negativeMap','pointPst','colorTmp', 'rawData_pn_meanOfROI','positiveRawAll','negativeRawAll')
figure(105)
saveas(gcf,'meanOfROI_positive','tiff')
figure(106)
saveas(gcf,'meanOfGroups_negative','tiff')
figure(107)
saveas(gcf,'meanOfGroups_positiveNegativeMix','tiff')
figure(109)
saveas(gcf,'distribution','tiff')
figure(108)
saveas(gcf,'PN','tiff')
figure(110)
saveas(gcf,'meanOfROI_abs','tiff')
imwrite(pnMap,'pnMap_raw.tif')
imwrite(positiveMap,'positiveMap.tif');
imwrite(negativeMap,'negativeMap.tif');
a=or(positiveMap,negativeMap);
imwrite(a,'pos_or_neg.tif')













% --- Executes on button press in pushbutton150.
function pushbutton150_Callback(hObject, eventdata, handles)
file=handles.file;
fileType=handles.fileType;
fileName=handles.fileName;
filePath=handles.filePath;
pointPst=handles.pointPst;
p=length(fileName);
cd(filePath)
img=differentTypeRead(file,fileType);
[m,n,nn_trash]=size(img);
groupN=length(pointPst);
ROI_N=zeros(groupN,1);
for ii=1:groupN
    ROI_N(ii)=length(pointPst{ii});
end
rawData_positive=zeros(p,sum(ROI_N(:))+1);
rawData_negative=zeros(p,sum(ROI_N(:))+1);
rawData_pn=zeros(p,sum(ROI_N(:))+1);
timeCourseTmp=eval(get(handles.et_timeCourse,'string'));
timeCourse=linspace(timeCourseTmp(1),timeCourseTmp(end),p).';
handles.timeCourse=timeCourse;
handles.timeCourseUnit=get(handles.et_unit,'string');
timeCourseUnit=handles.timeCourseUnit;
rawData_positive(:,1)=timeCourse;
rawData_negative(:,1)=timeCourse;
rawData_pn(:,1)=timeCourse;
positiveMap=handles.positiveMap;
negativeMap=handles.negativeMap;
pnYellow=and(positiveMap,handles.negativeMap);
positiveMap2=xor(pnYellow,positiveMap);
negativeMap2=xor(pnYellow,negativeMap);
positiveRawAll=zeros(p,sum(positiveMap2(:)),'single');
negativeRawAll=zeros(p,sum(negativeMap2(:)),'single');

h_wait=waitbar(0,'wait');
for kk=1:p
    waitbar(kk/p,h_wait,[num2str(100*kk/p,'%04.1f'),'%completed']);
    file=fileName(kk).name;
    img=differentTypeRead(file,fileType);
    positiveRawAll(kk,:)=(img(positiveMap2)).';
     negativeRawAll(kk,:)=(img(negativeMap2)).';
    if get(handles.cb_pnIOS_filter,'value')
        img=differentTypeReadFilter_handles(img,handles);
    end
    
    vv=1;
    for ii=1:length(pointPst)
        for jj=1:length(pointPst{ii})
            vv=vv+1;
            cx=pointPst{ii}{jj}(:,1);
            cy=pointPst{ii}{jj}(:,2);
            bwTmp=roipoly(img(:,:,1),cx,cy);
            bw_positive=and(bwTmp,positiveMap2);
            bw_negative=and(bwTmp,negativeMap2);
%             bw_positive=and(bwTmp,handles.positiveMap);
%             bw_negative=and(bwTmp,handles.negativeMap);
            bw_pn=and(bwTmp,pnYellow);
            ROItmp_positive=img(bw_positive);
            ROItmp_negative=img(bw_negative);
            ROItmp_pn=img(bw_pn);
            rawData_positive(kk,vv)=mean(ROItmp_positive(:));
            rawData_negative(kk,vv)=mean(ROItmp_negative(:));
            rawData_pn(kk,vv)=mean(ROItmp_pn(:));
        end
    end
end
close(h_wait)
rawData_positive(:,2:end)= ...
    IOS_time_gui_filter(rawData_positive(:,2:end),handles,'vertical');
rawData_negative(:,2:end)= ...
    IOS_time_gui_filter(rawData_negative(:,2:end),handles,'vertical');
rawData_pn(:,2:end)= ...
    IOS_time_gui_filter(rawData_pn(:,2:end),handles,'vertical');
handles.rawData_positive_meanOfROI=rawData_positive;
handles.rawData_negative_meanOfROI=rawData_negative;
handles.rawData_pn_meanOfROI=rawData_pn;
    handles.positiveRawAll=positiveRawAll;
     handles.negativeRawAll=negativeRawAll;
guidata(hObject,handles)
% hObject    handle to pushbutton150 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cb_pnIOS_filter.
function checkbox33_Callback(hObject, eventdata, handles)
% hObject    handle to cb_pnIOS_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_pnIOS_filter


% --- Executes on button press in pb_zstack_oneButton.

function pb_zstack_oneButton_Callback(hObject, eventdata, handles)
set(handles.pb_zstack_oneButton,'backgroundcolor',handles.bgColor);
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
% img=differentTypeRead(fullfile(filePath,file),fileType);
% ref=fft2(img);
planeNo=eval(get(handles.et_planeNo,'string'));
resultPathtrash=cell(planeNo,1);
fileNametrash=cell(planeNo,1);

for jj=1:planeNo
    result=['fftRegis_p',num2str(jj)];
    resultPath=fullfile(filePath,result);resultPathtrash{jj}=resultPath;
    fileNametrash{jj}=fileName(jj:planeNo:end);
    if exist(resultPath)==7
    else
        mkdir(resultPath)
    end
    for ii=jj:planeNo:p
        file=fileName(ii).name;
        img=differentTypeRead(fullfile(filePath,file),fileType);
        differentTypeWrite(img,fullfile(resultPath,file),fileType,handles.imgDepth);
    end
%     handles.filePath=resultPath;
%     cd(resultPath)
%     handles.file=fileName(ii).name;
%     
%     img=differentTypeRead(fullfile(filePath,file),fileType);
%     img=differentTypeReadFilter_handles(img,handles);
%     [output, img2] = dftregistrationIOS(ref,fft2(double(img)),a);
%     shift(ii,:)=output(3:4);
%     differentTypeWrite(abs(ifft2(img2)),fullfile(resultPath,file),fileType,handles.imgDepth);
end

handles.file01=file;
handles.fileName01=fileName;
handles.filePath01=filePath;
for jj=1:planeNo
    handles.fileName=fileNametrash{jj};
    handles.filePath=resultPathtrash{jj};
    handles=handlesUpdata(handles,eventdata);
    guidata(hObject,handles);
    %%
    handles=pushbutton135_Callback(hObject, eventdata, handles);
    handles=pushbutton153_Callback(hObject, eventdata, handles);
end
handles.fileName=handles.fileName01;
handles.filePath=handles.filePath01;
handles=handlesUpdata(handles,eventdata);
set(handles.pb_zstack_oneButton,'backgroundcolor',[0.678,0.922,1]);
guidata(hObject,handles)


% hObject    handle to pb_zstack_oneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function et_planeNo_Callback(hObject, eventdata, handles)
% hObject    handle to et_planeNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_planeNo as text
%        str2double(get(hObject,'String')) returns contents of et_planeNo as a double


% --- Executes during object creation, after setting all properties.
function et_planeNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_planeNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton153.
function handles=pushbutton153_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),fileType);
% ref=fft2(img);
result='Combined_tif_images';
result1 = 'Avg_mat_images';
cd(filePath);

filePath1=cd;

resultPath1=fullfile(fileparts(fileparts(fileparts(filePath1))),result);
mkdir(resultPath1)

resultPath2=fullfile(fileparts(fileparts(fileparts(filePath1))),result1);
mkdir(resultPath2)

resultPath=fullfile(fileparts(fileparts(filePath1)),result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end
filePath1=fileparts(filePath);

% result1='grayBackgroundImages';
% resultPath1=fullfile(filePath1,result1);
% mkdir(resultPath1);

cd(filePath1);


fileName1=dir('*.tif');

% aa=differentTypeRead(fullfile(filePath1,file),fileType);
aa=imread(fileName1(1).name);%fileName(1).name = a00000.BMP;
aa=aa(:,:,1);
[m,n]=size(aa);
img2=uint8(zeros(m,n));
imgTemp=double(zeros(m,n));

p2=size(fileName1,1);
loop=floor(p2/p);

%
h_wait=waitbar(0,'please wait');


D=uint8(zeros(m,n,3));
B=uint8(zeros(m,n,3));
AvgImg = single(zeros(m,n));
ImgTT=uint8(zeros(m,n,3));
D(:,:,1)=aa;
D(:,:,2)=aa;
D(:,:,3)=aa;


for tt=1:p
    %
    waitbar(tt/p,h_wait,[num2str(tt/p*100,'%04.1f'),'%completed']);
    for ii=((tt-1)*loop+1):(tt*loop)
        
        file1=fileName1(ii).name;
        img1=differentTypeRead(fullfile(filePath1,file1),'tif');
%         img1=differentTypeRead(fullfile(filePath1,file1),'mat');
        %     img=differentTypeReadFilter_handles(img,handles);
        imgTemp=double(imgTemp+img1);
    end
    img2=im2uint8(mat2gray(img1));
    
    file=fileName(tt).name;
    img=differentTypeRead(fullfile(filePath,file),fileType);
    AvgImg = AvgImg + img/p;
    
    newIntensity=inline(get(handles.et_inline,'string'));
    contents = cellstr(get(handles.pp_IOScolormap,'String'));
    colorSelectedTmp=contents(get(handles.pp_IOScolormap,'value'));
    colorSelected=eval(colorSelectedTmp{1});
    img3=newIntensity(img);
    imgColor=ind2rgb(uint8(img3),colorSelected);
    
    
    
    B=im2uint8(imgColor);
    
    D(:,:,1)=img2;
    D(:,:,2)=img2;
    D(:,:,3)=img2;
    for ss=1:m
        for zz=1:n
            if (B(ss,zz,3)==0 && B(ss,zz,2)<255) %%red only
            %if (B(ss,zz,1)==0 && B(ss,zz,2)<255)||(B(ss,zz,3)==0 && B(ss,zz,2)<255)  %%%blue||red
                ImgTT(ss,zz,:)=B(ss,zz,:);
            else
                ImgTT(ss,zz,:)=D(ss,zz,:);
                
            end
        end
    end
    
    cd(resultPath);
    indexNO=strfind(file1,'.tif');
    
    scaleName=get(handles.et_inline,'string');
    dotNO1=strfind(scaleName,'/');
    scale=scaleName(dotNO1+1:end);
%     imwrite(ImgTT,['combined_',file1(1:indexNO(end)-1),'.tif']);
 imwrite(ImgTT,['combined_scale',scale,'_',file1(1:indexNO(1)-1),'.tif']);
    cd(resultPath1)
 imwrite(ImgTT,['combined_scale',scale,'_',file1(1:indexNO(1)-1),'.tif']);
  cd(resultPath);
end

cd(resultPath2)
saveName=[file(1:end-3),'mat'];
save(saveName,'AvgImg');


handles.filePath=resultPath;
handles.fileType='tif';
cd(resultPath);
fileName=dir(['*.','tif']);
handles.fileName=fileName;
handles.file=['combined_scale',scale,'_',file1(1:indexNO(1)-1),'.tif'];
file=['combined_scale',scale,'_',file1(1:indexNO(1)-1),'.tif'];
handles.p=length(handles.fileName);
p=handles.p;

listName=cell(p,1);

for ii=1:p
    listName{ii}=handles.fileName(ii).name;
end
set(handles.dt_IntestedImg,'string',['1:',num2str(p)])
set(handles.lb_curFileNames,'string',listName)
set(handles.lb_curFileNames,'value',1)
handles=imageShow(handles);
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles)
 


function edit116_Callback(hObject, eventdata, handles)
% hObject    handle to edit116 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit116 as text
%        str2double(get(hObject,'String')) returns contents of edit116 as a double


% --- Executes during object creation, after setting all properties.
function edit116_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit116 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_zprj.
function pb_zprj_Callback(hObject, eventdata, handles)

set(hObject, 'BackgroundColor', 'r');
set(hObject, 'String','Projecting...');
drawnow;

file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;

file=fileName(1).name;
img=differentTypeRead(fullfile(filePath,file),fileType);
[m,n]=size(img);
img=zeros(m,n);
% img=differentTypeRead(fullfile(filePath,file),fileType);
% ref=fft2(img);
planeNo=eval(get(handles.et_planeNo,'string'));

result=['Bin_',num2str(planeNo),'mean'];
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end
% fileNametrash=cell(planeNo,1);
resultPathtrash=resultPath;
fileNametrash=cell(p/planeNo,1);

hw = waitbar(0,'wait','name','Z projection in progress');
for jj=1:p/planeNo
    waitbar(jj*planeNo/p,hw,[num2str(jj*planeNo/p*100,'%04.1f'),'%completed']);
    img=zeros(m,n);
%     fileNametrash{jj}=fileName(jj:planeNo:end);
    for ii=(1+(jj-1)*planeNo):(jj*planeNo)
        file=fileName(ii).name;
        img=img+differentTypeRead(fullfile(filePath,file),fileType);
    end
    img=img/planeNo;
    differentTypeWrite(img,fullfile(resultPath,file),fileType,handles.imgDepth);
end
close(hw);
handles.file01=file;
handles.fileName01=fileName;
handles.filePath01=filePath;

handles.filePath=resultPath;

% handles=pb_rawIOS_mat_Callback(hObject, eventdata, handles);
% handles.filePath=handles.mat_resultPath;
% handles.fileType='mat';
handles.fileName=dir(fullfile(handles.filePath,'*.tif'));
handles.file=handles.fileName(1).name;
handles.p=length(handles.fileName);
p=handles.p;
listName=cell(p,1);
for ii=1:p
    listName{ii}=handles.fileName(ii).name;
    
end
set(handles.dt_IntestedImg,'string',['1:',num2str(p)]);
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
handles=imageShow(handles);
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
set(hObject, 'BackgroundColor', [0.94,0.94,0.94]);
set(hObject, 'String','z_projection');
guidata(hObject,handles)
% guidata(hObject,handles);

%     handles.fileName=fileNametrash{jj};

% handles=handlesUpdata(handles,eventdata);
% guidata(hObject,handles);
%%
% handles=pushbutton135_Callback(hObject, eventdata, handles);
% handles=pushbutton153_Callback(hObject, eventdata, handles);
% 
% handles.fileName=handles.fileName01;
% handles.filePath=handles.filePath01;
% handles=handlesUpdata(handles,eventdata);
% guidata(hObject,handles)
% hObject    handle to pb_zprj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function et_shift_Callback(hObject, eventdata, handles)
% hObject    handle to et_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pointPst=handles.pointPst;
a=eval(get(handles.et_shift,'string'));
y1=a(1);
x1=a(2);
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
        cx=pointPst{ii}{jj}(:,1);
        handles.pointPst{ii}{jj}(:,1)=cx+x1;
        cy=pointPst{ii}{jj}(:,2);
        handles.pointPst{ii}{jj}(:,2)=cy+y1;

    end
end
pointPst=handles.pointPst;
pointPstH=cell(size(pointPst));
pointCurrent_h2=handles.pointCurrent_h;
handles.pointCurrent_h=[];
handles.pointCurrent_group=[];
ss=1;
for ii=1:length(pointPst)
    for jj=1:length(pointPst{ii})
%          delete(pointCurrent_h2(ss));
        if ishandle(pointCurrent_h2(ss))
            delete(pointCurrent_h2(ss));
        end
        ss=ss+1;
        
        cx=pointPst{ii}{jj}(:,1);
        cy=pointPst{ii}{jj}(:,2);
        hh=plot([cx;cx(1)],[cy;cy(1)],'color',handles.colorTmp(ii,:));
        handles.pointCurrent_h=[handles.pointCurrent_h;hh];
        handles.pointCurrent_group=[handles.pointCurrent_group;ii];
        pointPstH{ii}=[pointPstH{ii},hh];
    end
end
handles.pointPst=pointPst;
handles.pointPstH=pointPstH;
handles=GroupROITag(handles);
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of et_shift as text
%        str2double(get(hObject,'String')) returns contents of et_shift as a double


% --- Executes during object creation, after setting all properties.
function et_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_Your_background.
function pb_Your_background_Callback(hObject, eventdata, handles)
file=handles.fileMat;
filePath=handles.filePathMat;
fileName=handles.fileNameMat;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),'mat');
% ref=fft2(img);
filePath1=cd;

% result1='grayBackgroundImages';
% resultPath1=fullfile(filePath1,result1);
% mkdir(resultPath1);

cd(filePath);
cd ..

fileName1=dir('*.tif');



% aa=differentTypeRead(fullfile(filePath1,file),fileType);
aa=imread(fileName1(1).name);%fileName(1).name = a00000.BMP;
aa=aa(:,:,1);
[m,n]=size(aa);
img2=uint8(zeros(m,n));
imgTemp=double(zeros(m,n));

p2=size(fileName1,1);
loop=floor(p2/p);

%
h_wait=waitbar(0,'please wait');
%


D=uint8(zeros(m,n,3));
B=uint8(zeros(m,n,3));
ImgTT=uint8(zeros(m,n,3));
D(:,:,1)=aa;
D(:,:,2)=aa;
D(:,:,3)=aa;

[file00, filePath00]=uigetfile('*.*','please choose gray iamges');
img1=imread([filePath00,file00]);
img2=im2uint8(mat2gray(img1));

result='Combined_images_YourBGD';
resultPath=fullfile(handles.filePath_output,result);
if exist(resultPath)~=7
    mkdir(resultPath)
end


for tt=1:p
    waitbar(tt/p,h_wait,[num2str(tt/p*100,'%04.1f'),'%completed']);
%     for ii=((tt-1)*loop+1):(tt*loop)
%     for ii=1:loop
%         file1=fileName1(ii).name;
%         img1=differentTypeRead(fullfile(filePath1,file1),'tif');
% %         img1=differentTypeRead(fullfile(filePath1,file1),'mat');
%         %     img=differentTypeReadFilter_handles(img,handles);
%         imgTemp=double(imgTemp+img1);
%     end
    
   
    
    file=fileName(tt).name;
    img=differentTypeRead(fullfile(filePath,file),'mat');
    newIntensity=inline(get(handles.et_inline,'string'));
    contents = cellstr(get(handles.pp_IOScolormap,'String'));
    colorSelectedTmp=contents(get(handles.pp_IOScolormap,'value'));
    colorSelected=eval(colorSelectedTmp{1});
    img3=newIntensity(img);
    imgColor=ind2rgb(uint8(img3),colorSelected);
    
    
    
    B=im2uint8(imgColor);
 
    D(:,:,1)=img2;
    D(:,:,2)=img2;
    D(:,:,3)=img2;
    
    for ss=1:m
        for zz=1:n
            if (B(ss,zz,1)==0 && B(ss,zz,2)<255)||(B(ss,zz,3)==0 && B(ss,zz,2)<255)  %%%blue||red
                ImgTT(ss,zz,:)=B(ss,zz,:);
            else
                ImgTT(ss,zz,:)=D(ss,zz,:);
                
            end
        end
    end
    
    cd(resultPath);
    indexNO=strfind(file00,'.tif');
    
    scaleName=get(handles.et_inline,'string');
    dotNO1=strfind(scaleName,'/');
    scale=scaleName(dotNO1+1:end);
%     imwrite(ImgTT,['combined_',file1(1:indexNO(end)-1),'.tif']);
 imwrite(ImgTT,['yb_scale',scale,'_',file00(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif']);
end

handles.filePath=resultPath;
handles.fileType='tif';
cd(resultPath);
fileName=dir(['*.','tif']);
handles.fileName=fileName;
handles.file=['yb_scale',scale,'_',file00(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif'];
file=['yb_scale',scale,'_',file00(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif'];
handles.p=length(handles.fileName);
p=handles.p;

listName=cell(p,1);

for ii=1:p
    listName{ii}=handles.fileName(ii).name;
end
set(handles.dt_IntestedImg,'string',['1:',num2str(p)]);
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
handles=imageShow(handles);
close(h_wait);
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles)
% hObject    handle to pb_Your_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_Your_background_NOBlue.
function pb_Your_background_NOBlue_Callback(hObject, eventdata, handles)
file=handles.fileMat;
filePath=handles.filePathMat;
fileName=handles.fileNameMat;
p=length(fileName);
fileType=handles.fileType;
img=differentTypeRead(fullfile(filePath,file),'mat');
% ref=fft2(img);

filePath1=cd;

% result1='grayBackgroundImages';
% resultPath1=fullfile(filePath1,result1);
% mkdir(resultPath1);

cd(filePath);
cd ..

fileName1=dir('*.tif');


% aa=differentTypeRead(fullfile(filePath1,file),fileType);
aa=imread(fileName1(1).name);%fileName(1).name = a00000.BMP;
aa=aa(:,:,1);
[m,n]=size(aa);
img2=uint8(zeros(m,n));
imgTemp=double(zeros(m,n));

p2=size(fileName1,1);
loop=p2/p;

%
h_wait=waitbar(0,'please wait');
%


D=uint8(zeros(m,n,3));
B=uint8(zeros(m,n,3));
ImgTT=uint8(zeros(m,n,3));
D(:,:,1)=aa;
D(:,:,2)=aa;
D(:,:,3)=aa;

[file00, filePath00]=uigetfile('*.*','please choose gray iamges');
img1=imread([filePath00,file00]);
img2=im2uint8(mat2gray(img1));


result='Combined_images_YourBGD_noBlue';
resultPath=fullfile(handles.filePath_output,result);

if exist(resultPath)~=7
    mkdir(resultPath)
end

for tt=1:p
    waitbar(tt/p,h_wait,[num2str(tt/p*100,'%04.1f'),'%completed']);
    %     for ii=((tt-1)*loop+1):(tt*loop)
    %     for ii=1:loop
    %         file1=fileName1(ii).name;
    %         img1=differentTypeRead(fullfile(filePath1,file1),'tif');
    % %         img1=differentTypeRead(fullfile(filePath1,file1),'mat');
    %         %     img=differentTypeReadFilter_handles(img,handles);
    %         imgTemp=double(imgTemp+img1);
    %     end
    
    
    
      file=fileName(tt).name;
    img=differentTypeRead(fullfile(filePath,file),'mat');
    newIntensity=inline(get(handles.et_inline,'string'));
    contents = cellstr(get(handles.pp_IOScolormap,'String'));
    colorSelectedTmp=contents(get(handles.pp_IOScolormap,'value'));
    colorSelected=eval(colorSelectedTmp{1});
    img3=newIntensity(img);
    imgColor=ind2rgb(uint8(img3),colorSelected);
    
    
    B=im2uint8(imgColor);
    
    D(:,:,1)=img2;
    D(:,:,2)=img2;
    D(:,:,3)=img2;
    for ss=1:m
        for zz=1:n
            if (B(ss,zz,3)==0 && B(ss,zz,2)<255)  %%%blue||red
                ImgTT(ss,zz,:)=B(ss,zz,:);
            else
                ImgTT(ss,zz,:)=D(ss,zz,:);
                
            end
        end
    end
    
    
    
    cd(resultPath);
    indexNO=strfind(file00,'.tif');
    
    scaleName=get(handles.et_inline,'string');
    dotNO1=strfind(scaleName,'/');
    scale=scaleName(dotNO1+1:end);
    %     imwrite(ImgTT,['combined_',file1(1:indexNO(end)-1),'.tif']);
    imwrite(ImgTT,['ybnb_scale',scale,'_',file00(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif']);
end

handles.filePath=resultPath;
handles.fileType='tif';
cd(resultPath);
fileName=dir(['*.','tif']);
handles.fileName=fileName;
handles.file=['ybnb_scale',scale,'_',file00(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif'];
file=['ybnb_scale',scale,'_',file00(1:indexNO(1)-5),'_',num2str(tt,'%03d'),'.tif'];
handles.p=length(handles.fileName);
p=handles.p;

listName=cell(p,1);

for ii=1:p
    listName{ii}=handles.fileName(ii).name;
end
set(handles.dt_IntestedImg,'string',['1:',num2str(p)]);
set(handles.lb_curFileNames,'string',listName);
set(handles.lb_curFileNames,'value',1);
set(handles.st_curFilePath,'string',handles.filePath);
handles=imageShow(handles);
close(h_wait);
handles=lb_curFileNames_Callback(handles.lb_curFileNames, eventdata, handles);
guidata(hObject,handles)
% hObject    handle to pb_Your_background_NOBlue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cb_cdCurrent.
function cd_cdCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to cb_cdCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_cdCurrent


% --- Executes on button press in pushbutton157.
function pushbutton157_Callback(hObject, eventdata, handles)
file=handles.file;
filePath=handles.filePath;
fileName=handles.fileName;
p=length(fileName);
fileType=handles.fileType;
% img=differentTypeRead(fullfile(filePath,file),fileType);
% ref=fft2(img);
planeNo=eval(get(handles.et_planeNo,'string'));
resultPathtrash=cell(planeNo,1);
fileNametrash=cell(planeNo,1);

for jj=1:planeNo
    result=['fftRegis_p',num2str(jj)];
    resultPath=fullfile(filePath,result);resultPathtrash{jj}=resultPath;
    fileNametrash{jj}=fileName(jj:planeNo:end);
    if exist(resultPath)==7
    else
        mkdir(resultPath)
    end
    for ii=jj:planeNo:p
        file=fileName(ii).name;
        img=differentTypeRead(fullfile(filePath,file),fileType);
        differentTypeWrite(img,fullfile(resultPath,file),fileType,handles.imgDepth);
    end
%     handles.filePath=resultPath;
%     cd(resultPath)
%     handles.file=fileName(ii).name;
%     
%     img=differentTypeRead(fullfile(filePath,file),fileType);
%     img=differentTypeReadFilter_handles(img,handles);
%     [output, img2] = dftregistrationIOS(ref,fft2(double(img)),a);
%     shift(ii,:)=output(3:4);
%     differentTypeWrite(abs(ifft2(img2)),fullfile(resultPath,file),fileType,handles.imgDepth);
end

handles.file01=file;
handles.fileName01=fileName;
handles.filePath01=filePath;
for jj=1:planeNo
    handles.fileName=fileNametrash{jj};
    handles.filePath=resultPathtrash{jj};
    handles=handlesUpdata(handles,eventdata);
    guidata(hObject,handles);
    %%
    handles=pushbutton134_Callback(hObject, eventdata, handles);
    
end
handles.fileName=handles.fileName01;
handles.filePath=handles.filePath01;
handles=handlesUpdata(handles,eventdata);
guidata(hObject,handles)

% hObject    handle to pushbutton157 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over lb_curFileNames.
function lb_curFileNames_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to lb_curFileNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function imgQuality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgQuality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Data',cell(3,1));


% --- Executes during object creation, after setting all properties.
function axes_rawImg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_rawImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% appHome = fileparts(which(mfilename));
% % handles.axes_rawImg = hObject;
% % handles.filePath = [appHome,filesep];
% % handles.file = 'zfhc.png';
% imshow(imread([appHome,filesep,'zfhc.png']));
% % handles.fileType = 'png';
% hold off;
% axis off;
% guidata(hObject, handles)

% Hint: place code in OpeningFcn to populate axes_rawImg


% --- Executes on button press in pb_savePath.
function pb_savePath_Callback(hObject, eventdata, handles)
outPath = uigetdir;
handles.filePath_output = outPath;
set(handles.et_IOSresultPath,'String',outPath);
guidata(hObject,handles)

% hObject    handle to pb_savePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pb_oneButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pb_oneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on lb_curFileNames and none of its controls.
function lb_curFileNames_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to lb_curFileNames (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
