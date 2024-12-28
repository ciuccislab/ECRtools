function varargout = ECRTOOLS_GUI(varargin)
% ECRTOOLS_GUI MATLAB code for ECRTOOLS_GUI.fig
%      ECRTOOLS_GUI, by itself, creates a new ECRTOOLS_GUI or raises the existing
%      singleton*.
%
%      H = ECRTOOLS_GUI returns the handle to a new ECRTOOLS_GUI or the handle to
%      the existing singleton*.
%
%      ECRTOOLS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ECRTOOLS_GUI.M with the given input arguments.
%
%      ECRTOOLS_GUI('Property','Value',...) creates a new ECRTOOLS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ECRTOOLS_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ECRTOOLS_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help ECRTOOLS_GUI

% Last Modified by GUIDE v2.5 08-Nov-2012 17:53:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ECRTOOLS_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ECRTOOLS_GUI_OutputFcn, ...
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


% --- Executes just before ECRTOOLS_GUI is made visible.
function ECRTOOLS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ECRTOOLS_GUI (see VARARGIN)

% Choose default command line output for ECRTOOLS_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ECRTOOLS_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = ECRTOOLS_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function k_value_Callback(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of k_value as text
%        str2double(get(hObject,'String')) returns contents of k_value as a double

k_value = str2double(get(hObject, 'String'));
if isnan(k_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (k_value)<=0
    set(hObject, 'String', 0);
    errordlg('k needs to be a positive quantity','Error');
end

% Save the new k value
% note input is in cm/s 
% so need to multiply by 1E-2
handles.k_value = k_value*1E-2;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function k_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function D_value_Callback(hObject, eventdata, handles)
% hObject    handle to D_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D_value as text
%        str2double(get(hObject,'String')) returns contents of D_value as a double

D_value = str2double(get(hObject, 'String'));

if isnan(D_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (D_value)<=0
    set(hObject, 'String', 0);
    errordlg('D needs to be a positive quantity','Error');
end

% Save the new D value
% note input is in cm^2/s 
% so need to multiply by 1E-4
handles.D_value = D_value*1E-4;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function D_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when panel_k_D is resized.
function panel_k_D_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to panel_k_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function N_value_Callback(hObject, eventdata, handles)
% hObject    handle to N_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_value as text
%        str2double(get(hObject,'String')) returns contents of N_value as a double

N_value = str2double(get(hObject, 'String'))
if isnan(N_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (N_value)<=0
    set(hObject, 'String', 0);
    errordlg('Number of measurements needs to be a positive integer','Error');
elseif abs(floor(N_value)-N_value)>0
    set(hObject, 'String', 0);
    errordlg('Number of measurements needs to be an integer','Error');
end

% Save the new N value
handles.N_value = N_value;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function N_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_exp_value_Callback(hObject, eventdata, handles)
% hObject    handle to t_exp_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_exp_value as text
%        str2double(get(hObject,'String')) returns contents of t_exp_value as a double

t_exp_value = str2double(get(hObject, 'String'));
if isnan(t_exp_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (t_exp_value)<=0
    set(hObject, 'String', 0);
    errordlg('Time of experiment must be positive','Error');
elseif (t_exp_value)>5E5
    set(hObject, 'String', 0);
    errordlg('Time of experiment is greater than 5.5 days!','Error');
end

% Save the new t_exp value
handles.t_exp_value = t_exp_value;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function t_exp_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_exp_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function st_dev_value_Callback(hObject, eventdata, handles)
% hObject    handle to st_dev_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_dev_value as text
%        str2double(get(hObject,'String')) returns contents of st_dev_value as a double

st_dev_value = str2double(get(hObject, 'String'));
if isnan(st_dev_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (st_dev_value)<0
    set(hObject, 'String', 0);
    errordlg('Measurement error must be positive','Error');
elseif (st_dev_value)>50
    set(hObject, 'String', 0);
    errordlg('Measurement error is greater that 50%','Error');
elseif (st_dev_value)==0
    set(hObject, 'String', 0);
    errordlg('No error?','Error');
end

% Save the new standard deviation value
handles.st_dev_value = st_dev_value*1E-2;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function st_dev_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to st_dev_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lz_value_Callback(hObject, eventdata, handles)
% hObject    handle to Lz_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lz_value as text
%        str2double(get(hObject,'String')) returns contents of Lz_value as a double

Lz_value = str2double(get(hObject, 'String'));
if isnan(Lz_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (Lz_value)<0
    set(hObject, 'String', 0);
    errordlg('Sample Thickness error must be positive','Error');
elseif (Lz_value)>20 && Lz_value~=Inf
    set(hObject, 'String', 0);
    errordlg('The sample size is greater than 20 cm!','Error');
elseif (Lz_value)==0
    set(hObject, 'String', 0);
    errordlg('Does the sample exist?','Error');
end

% Save the new thickness value
% note input is in mm
% so need to multiply by 1E-2
% note also that ECR uses half thickness so you need to divide by two
handles.Lz_value = Lz_value*1E-2;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Lz_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lz_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when panel_thickness is resized.
function panel_thickness_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to panel_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in draw_ECR_push.
function draw_ECR_push_Callback(hObject, eventdata, handles)
% hObject    handle to draw_ECR_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.ECR_response);
cla;

a = ecr_base;

a.k_ref = handles.k_value; % value for k
a.D_ref = handles.D_value; % value for D

 
a.half_thickness_x_ref = 0.5*handles.Lx_value; % thickness in x direction
a.half_thickness_y_ref = 0.5*handles.Ly_value; % thickness in y direction
a.half_thickness_z_ref = 0.5*handles.Lz_value; % thickness in z_direction

a.standard_dev = handles.st_dev_value; % standard deviation
N_exp = handles.N_value;
t_exp = handles.t_exp_value; 

time = linspace(0, t_exp, N_exp);

sigma_n         = a.sigma_n_det(time);
sigma_n_meas    = a.sigma_n_meas(time);
error_meas      = sigma_n_meas - sigma_n;


[AX,H1,H2] = plotyy(time, sigma_n, time, 100*error_meas,'plot','plot');
set(H1,'LineStyle','-', 'Linewidth', 2, 'Color', 'k');
set(H2,'LineStyle','-.', 'MarkerSize', 12, 'Color', 'b');

set(AX(1),'ycolor','k')
set(AX(1),'FontSize',15)
set(AX(1),'XLim',[0 max(time)])
set(AX(1),'YLim',[0 1.0])
set(AX(1),'YTick',[0:0.1:1.0])
set(AX(2),'ycolor','b')
set(AX(2),'FontSize',15)
set(AX(2),'XLim',[0 max(time)])
set(AX(2),'YLim',[-100 100])
set(AX(2),'YTick',[-100:10:100])

hold on
plot(time, sigma_n_meas, '.r', 'MarkerSize', 15)
hold off
xlabel('$t_{\rm exp}/{\rm s}$', 'Interpreter' , 'Latex', 'Fontsize', 18)
%set(get(AX(1),'Ylabel'),'String','$\sigma_n$', 'Interpreter' , 'Latex', 'Fontsize', 18)
set(get(AX(1),'Ylabel'),'String','$\sigma_n$', 'Interpreter' , 'Latex', 'Fontsize', 18)
set(get(AX(2),'Ylabel'),'String','$100\times(\sigma_n-\sigma_n^{\rm meas}~)$', 'Interpreter' , 'Latex', 'Fontsize', 18)


% --- Executes on button press in draw_sensitivity_push.
function draw_sensitivity_push_Callback(hObject, eventdata, handles)
% hObject    handle to draw_sensitivity_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.ECR_sensitivity);
cla;

% define class
% use sensitivity class
a = ecr_sens;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Class Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

a.k_ref = handles.k_value; % value for k
a.D_ref = handles.D_value; % value for D

a.half_thickness_x_ref = 0.5*handles.Lx_value; % thickness in x direction
a.half_thickness_y_ref = 0.5*handles.Ly_value; % thickness in y direction
a.half_thickness_z_ref = 0.5*handles.Lz_value; % thickness in z_direction

% Define Timespan
N = 1000;
t_min = 1E-2; % Minimum time is 1E-2 s
t_max = 1E6; % Maximum time is 1E4 s
% time needs to be column vector
time = logspace(log10(t_min), log10(t_max), 1000)'; % Maximum time is 100 s

grad_out = a.grad_sigma_n(time);

semilogx(time, grad_out(:, 1)*100, '-k', 'LineWidth', 2)
hold on
semilogx(time, grad_out(:, 2)*100, '-r', 'LineWidth', 2)
max_grad = [grad_out(:, 1); grad_out(:, 2)];
max_grad = max(max_grad);

vec_dummy   = [handles.t_exp_value, handles.t_exp_value];
vec_dummy2  = vec_dummy/handles.N_value;
%semilogx(vec_dummy, 100*[0, max_grad], ':g', 'LineWidth', 2)
%semilogx(vec_dummy2, 100*[0, max_grad], ':g', 'LineWidth', 2)
semilogx(vec_dummy, [0, 50], ':g', 'LineWidth', 2)
semilogx(vec_dummy2, [0, 50], ':g', 'LineWidth', 2)

hold off
axis([t_min, t_max, 0, 40]);

xlabel('$t$/s', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$100\times q\frac{\partial \sigma_n}{\partial q}$', 'Interpreter', 'Latex', 'FontSize', 18)
set(gca,'FontSize',14)
legend('q=k', 'q=D', 'exp. range', 'Location', 'NorthWest')


% --- Executes on button press in draw_confidence_push.
function draw_confidence_push_Callback(hObject, eventdata, handles)
% hObject    handle to draw_confidence_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.ECR_confidence);
cla;

% define class
% use sensitivity class
a = ecr_sens;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Class Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

a.k_ref = handles.k_value; % value for k
a.D_ref = handles.D_value; % value for D

a.half_thickness_x_ref = 0.5*handles.Lx_value; % thickness in x direction
a.half_thickness_y_ref = 0.5*handles.Ly_value; % thickness in y direction
a.half_thickness_z_ref = 0.5*handles.Lz_value; % thickness in z_direction

a.standard_dev = handles.st_dev_value;

% Define Timespan
N_exp = handles.N_value;
t_exp = handles.t_exp_value; 
% time needs to be column vector
time = linspace(0, t_exp, N_exp)';

boundary_cov  = a.boundary_cov_sigma_n(time);

plot(boundary_cov(1,:), boundary_cov(2,:), '-k', 'LineWidth', 2)

xlabel('$\hat k/k_{\rm exact}$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$\hat D/D_{\rm exact}$', 'Interpreter', 'Latex', 'FontSize', 15)
set(gca,'FontSize',14)
axis('equal')
Delta_x = handles.size_confidence_value;

axis([1.-Delta_x, 1.+Delta_x, 1.-Delta_x, 1.+Delta_x])



function size_confidence_value_Callback(hObject, eventdata, handles)
% hObject    handle to size_confidence_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_confidence_value as text
%        str2double(get(hObject,'String')) returns contents of size_confidence_value as a double

size_confidence_value = str2double(get(hObject, 'String'));
if isnan(size_confidence_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (size_confidence_value)<=0
    set(hObject, 'String', 0);
    errordlg('Measurement error must be positive','Error');
end

% Save the new standard deviation value
handles.size_confidence_value = size_confidence_value*1E-2;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function size_confidence_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_confidence_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_data_push.
function load_data_push_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.fileName = uigetfile('*');
guidata(hObject, handles);


% --- Executes on button press in fit_data_push.
function fit_data_push_Callback(hObject, eventdata, handles)
% hObject    handle to fit_data_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fileName = handles.fileName;
data = load(fileName);
time_exp = data(:,1)-data(1,1);
sigma_n_exp = data(:,2);

a = ecr_fit;
a.k_ref = 1E-4; % k - unit of k is m/s
a.D_ref = 1E-9; % D - unit of D m^2/s

a.half_thickness_x_ref = 0.5*handles.Lx_value; % half thickness in x direction
a.half_thickness_y_ref = 0.5*handles.Ly_value; % half thickness in y direction
a.half_thickness_z_ref = 0.5*handles.Lz_value; % half thickness in z_direction

a.meas_time = time_exp;
a.meas_sigma_n = sigma_n_exp;

% first run fit within a large span
theta_out = a.log10_fit();

% update the reference value
a.k_ref = theta_out(1)*a.k_ref;
a.D_ref = theta_out(2)*a.D_ref;

theta_out = a.fit(handles.solver);
a.k_ref = theta_out(1)*a.k_ref;
a.D_ref = theta_out(2)*a.D_ref;

sigma_computed = a.output_sigma_n_det(theta_out);
error_meas = (sigma_computed-sigma_n_exp);
std_computed = a.compute_std_meas();
a.standard_dev = std_computed;

handles.k_value = a.k_ref;
k_string = num2str(1E2*a.k_ref);
myHandle = findobj('Tag','k_value');
set(myHandle, 'string', k_string);

handles.D_value = a.D_ref;
D_string = num2str(1E4*a.D_ref);
myHandle = findobj('Tag','D_value');
set(myHandle, 'string', D_string);

handles.N_value = numel(time_exp);
N_string = num2str(numel(time_exp));
myHandle = findobj('Tag','N_value');
set(myHandle, 'string', N_string);

handles.t_exp_value = max(time_exp);
t_exp_string = num2str(max(time_exp));
myHandle = findobj('Tag','t_exp_value');
set(myHandle, 'string', t_exp_string);

handles.st_dev_value = a.standard_dev;
st_dev_string = num2str(100*a.standard_dev);
myHandle = findobj('Tag','st_dev_value');
set(myHandle, 'string', st_dev_string);


% push the new handles 
% into the guidata 
guidata(hObject,handles); 
% otherwise handles is only
% a local variable

axes(handles.ECR_response);
cla;

[AX,H1,H2] = plotyy(time_exp, sigma_computed, time_exp, 100*error_meas,'plot','plot');
set(H1,'LineStyle','-', 'Linewidth', 2, 'Color', 'k');
set(H2,'Marker','.', 'MarkerSize', 12, 'Color', 'b');

set(AX(1),'ycolor','k')
set(AX(1),'FontSize',15)
set(AX(1),'XLim',[0 max(time_exp)])
set(AX(1),'YLim',[0 1.0])
set(AX(1),'YTick',[0:0.1:1.0])
set(AX(2),'ycolor','b')
set(AX(2),'FontSize',15)
set(AX(2),'XLim',[0 max(time_exp)])
set(AX(2),'YLim',[-100 100])
set(AX(2),'YTick',[-100:10:100])

hold on
plot(time_exp, sigma_n_exp, '.r', 'MarkerSize', 15)
hold off
xlabel('$t_{\rm exp}/{\rm s}$', 'Interpreter' , 'Latex', 'Fontsize', 18)
%set(get(AX(1),'Ylabel'),'String','$\sigma_n$', 'Interpreter' , 'Latex', 'Fontsize', 18)
set(get(AX(1),'Ylabel'),'String','$\sigma_n$', 'Interpreter' , 'Latex', 'Fontsize', 18)
set(get(AX(2),'Ylabel'),'String','$100\times(\sigma_n-\sigma_n^{\rm meas}~)$', 'Interpreter' , 'Latex', 'Fontsize', 18)





% --- Executes on selection change in solver_list.
function solver_list_Callback(hObject, eventdata, handles)
% hObject    handle to solver_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns solver_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from solver_list


% --- Executes during object creation, after setting all properties.
function solver_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solver_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_state_gui_push.
function save_state_gui_push_Callback(hObject, eventdata, handles)
% hObject    handle to save_state_gui_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveState(handles);

% --- Executes on button press in load_state_gui_push.
function load_state_gui_push_Callback(hObject, eventdata, handles)
% hObject    handle to load_state_gui_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

restoreState(hObject, handles);


function saveState(handles)

state.k_value               = handles.k_value;
state.D_value               = handles.D_value;
state.N_value               = handles.N_value;
state.t_exp_value           = handles.t_exp_value;
state.Lx_value              = handles.Lx_value;
state.Ly_value              = handles.Ly_value;
state.Lz_value              = handles.Lz_value;
state.st_dev_value          = handles.st_dev_value;
state.size_confidence_value = handles.size_confidence_value;

save state.mat state

function restoreState(hObject, handles)

fileName='state.mat';

if exist(fileName)
    
    load(fileName);

    handles.k_value = state.k_value;
    k_string = num2str(1E2*state.k_value);
    myHandle = findobj('Tag','k_value');
    set(myHandle, 'string', k_string);

    handles.D_value = state.D_value;
    D_string = num2str(1E4*state.D_value);
    myHandle = findobj('Tag','D_value');
    set(myHandle, 'string', D_string);

    handles.N_value = state.N_value;
    N_string = num2str(state.N_value);
    myHandle = findobj('Tag','N_value');
    set(myHandle, 'string', N_string);

    handles.t_exp_value = state.t_exp_value;
    t_exp_string = num2str(state.t_exp_value);
    myHandle = findobj('Tag','t_exp_value');
    set(myHandle, 'string', t_exp_string);

    handles.Lx_value = state.Lx_value;
    Lx_string = num2str(1E2*state.Lx_value);
    myHandle = findobj('Tag','Lx_value');
    set(myHandle, 'string', Lx_string);
    
    handles.Ly_value = state.Ly_value;
    Ly_string = num2str(1E2*state.Ly_value);
    myHandle = findobj('Tag','Ly_value');
    set(myHandle, 'string', Ly_string);

    handles.Lz_value = state.Lz_value;
    Lz_string = num2str(1E2*state.Lz_value);
    myHandle = findobj('Tag','Lz_value');
    set(myHandle, 'string', Lz_string);
    
    handles.st_dev_value = state.st_dev_value;
    st_dev_string = num2str(100*state.st_dev_value);
    myHandle = findobj('Tag','st_dev_value');
    set(myHandle, 'string', st_dev_string);

    handles.size_confidence_value = state.size_confidence_value;
    size_confidence_string = num2str(100*state.size_confidence_value);
    myHandle = findobj('Tag','size_confidence_value');
    set(myHandle, 'string', size_confidence_string);
    
end

% push the new handles 
% into the guidata 
guidata(hObject,handles); 
% otherwise handles in only
% a local variable


% --- Executes on button press in populate_confidence_push.
function populate_confidence_push_Callback(hObject, eventdata, handles)
% hObject    handle to populate_confidence_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ECR_confidence);
cla;

% Define Timespan
N_exp = handles.N_value;
t_exp = handles.t_exp_value; 
% time needs to be column vector
time = linspace(0, t_exp, N_exp)';

% use sensitivity class
a = ecr_fit;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Class Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

a.k_ref = handles.k_value; % value for k
a.D_ref = handles.D_value; % value for D

a.half_thickness_x_ref = 0.5*handles.Lx_value; % thickness in x direction
a.half_thickness_y_ref = 0.5*handles.Ly_value; % thickness in y direction
a.half_thickness_z_ref = 0.5*handles.Lz_value; % thickness in z_direction

a.standard_dev = handles.st_dev_value;
a.meas_time = time; % experimental time - unit is s

boundary_cov  = a.boundary_cov_sigma_n(time);
theta_vec = a.synthetic_exp(handles.solver, handles.N_points_confidence_value);

plot(boundary_cov(1,:), boundary_cov(2,:), '-k', 'LineWidth', 2)
hold on
plot(theta_vec(1,:), theta_vec(2,:), '.r', 'LineWidth', 2)

xlabel('$\hat k/k_{\rm exact}$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$\hat D/D_{\rm exact}$', 'Interpreter', 'Latex', 'FontSize', 15)
set(gca,'FontSize',14)
axis('equal')
Delta_x = handles.size_confidence_value;

axis([1.-Delta_x, 1.+Delta_x, 1.-Delta_x, 1.+Delta_x]);


% --- Executes on selection change in select_solver.
function select_solver_Callback(hObject, eventdata, handles)
% hObject    handle to select_solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_solver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_solver
% 
% can do either 
% contents = cellstr(get(hObject,'String'))
% contents{get(hObject,'Value')}
% 
% or conversely
index_selected = get(hObject,'Value');
list = get(hObject,'String');
item_selected = list{index_selected}; % Convert from cell array to string

handles.solver = item_selected;
% push the new handles 
% into the guidata 
guidata(hObject,handles); 
% otherwise handles in only
% a local variable



% --- Executes during object creation, after setting all properties.
function select_solver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_points_confidence_value_Callback(hObject, eventdata, handles)
% hObject    handle to N_points_confidence_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_points_confidence_value as text
%        str2double(get(hObject,'String')) returns contents of N_points_confidence_value as a double

N_points_confidence_value = str2double(get(hObject, 'String'))
if isnan(N_points_confidence_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (N_points_confidence_value)<=0
    set(hObject, 'String', 0);
    errordlg('Number of measurements needs to be a positive integer','Error');
elseif abs(floor(N_points_confidence_value)-N_points_confidence_value)>0
    set(hObject, 'String', 0);
    errordlg('Number of measurements needs to be an integer','Error');
end

% Save the new N value
handles.N_points_confidence_value = N_points_confidence_value;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function N_points_confidence_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_points_confidence_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stoch_points_txt_Callback(hObject, eventdata, handles)
% hObject    handle to stoch_points_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stoch_points_txt as text
%        str2double(get(hObject,'String')) returns contents of stoch_points_txt as a double


% --- Executes during object creation, after setting all properties.
function stoch_points_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stoch_points_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lx_value_Callback(hObject, eventdata, handles)
% hObject    handle to Lx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lx_value as text
%        str2double(get(hObject,'String')) returns contents of Lx_value as a double
% hObject    handle to Lz_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lz_value as text
%        str2double(get(hObject,'String')) returns contents of Lz_value as a double

Lx_value = str2double(get(hObject, 'String'));
if isnan(Lx_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (Lx_value)<0
    set(hObject, 'String', 0);
    errordlg('Sample Thickness error must be positive','Error');
elseif (Lx_value)>20
    set(hObject, 'String', 0);
    errordlg('The sample size is greater than 20 cm!','Error');
elseif (Lx_value)==0
    set(hObject, 'String', 0);
    errordlg('Does the sample exist?','Error');
end

% Save the new thickness value
% note input is in mm
% so need to multiply by 1E-2
% note also that ECR uses half thickness so you need to divide by two
handles.Lx_value = Lx_value*1E-2;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Lx_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ly_value_Callback(hObject, eventdata, handles)
% hObject    handle to Ly_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ly_value as text
%        str2double(get(hObject,'String')) returns contents of Ly_value as a double
% hObject    handle to Lz_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lz_value as text
%        str2double(get(hObject,'String')) returns contents of Lz_value as a double

Ly_value = str2double(get(hObject, 'String'));
if isnan(Ly_value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif (Ly_value)<0
    set(hObject, 'String', 0);
    errordlg('Sample Thickness error must be positive','Error');
elseif (Ly_value)>20 && Ly_value~=Inf
    set(hObject, 'String', 0);
    errordlg('The sample size is greater than 20 cm!','Error');
elseif (Ly_value)==0
    set(hObject, 'String', 0);
    errordlg('Does the sample exist?','Error');
end

% Save the new thickness value
% note input is in mm
% so need to multiply by 1E-2
% note also that ECR uses half thickness so you need to divide by two
handles.Ly_value = Ly_value*1E-2;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Ly_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ly_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
