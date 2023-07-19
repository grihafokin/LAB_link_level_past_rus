% Фокин Г.А., СПбГУТ, 2023.
% имитационная модель оценки помех для двух радиолиний gNB-UE
% с диаграммообразованием на основе позиционирования 
% при использовании на базовых станциях прямоугольных антенных решеток (АР)
% имитационная модель использует встроенные функции пакета расширения 
% Phased Array System Toolbox (PAST)
% For citation: Fokin G. Location aware beamforming in mm-wave band 
% ultra-dense radio access networks. Part 1. Model of two links. 
% Proc. of Telecom. Universities. 2023 XX (in Russ.) DOI:XXX
clear all; close all; clc; 
% set(0,'DefaultFigureWindowStyle','docked');
plot_3D = 1;      % построение графиков в 3D
record_video = 1; % запись анимации работы модели; установить plot_3D = 1

%% инициализация параметров прямоугольной антенной решетки на gNB
fc = 30e9;                   % несущая частота, Гц
c = physconst('LightSpeed'); % скорость света, м/с
lambda = c/fc;               % длина волны, м 
drow = lambda/2;             % расстояние между строками АР, м
dcol = lambda/2;             % расстояние между столбцами АР, м 
nrow = 8;                    % число строк в массиве АР, штук 
ncol = 8;                    % число столбцов в массиве АР, штук 
% формирование прямоугольной АР 8×8 URA1 на gNB1
URA1 = phased.URA(...
    'Size',           [nrow ncol], ...
    'ElementSpacing', [drow dcol], ...
    'ArrayNormal',    'y');
URA1.Element.BackBaffled = true;
% формирование прямоугольной АР 8×8 URA2 на gNB2
URA2 = phased.URA(...
    'Size',           [nrow ncol], ...
    'ElementSpacing', [drow dcol], ...
    'ArrayNormal',    'y');
URA2.Element.BackBaffled = true;
% оценка визуализация параметров прямоугольной АР
% beamwidth(URA1,fc,'dBDown',3); % вычисление/построение ширины луча АР
% pattern(URA1,fc);            % построение ДНА
% patternAzimuth(URA1,fc);     % построение ДНА в горизонтальной плоскости
% patternElevation(URA1,fc);   % построение ДНА в вертикальной плоскости
% viewArray(URA1);view(45,45); % визуализация массива элементов АР

%% инициализация параметров территориального распределения gNB и UE
traffic_light_h = [0; 0; 0]; % коэффициенты 3D STL объекта светофора
% координаты URA и начальные азимуты и углы места 
URA1_pos = [10; 0; 8] + traffic_light_h; % координаты АР URA1 на gNB1
URA2_pos = [20; 0; 8] + traffic_light_h; % координаты АР URA2 на gNB2
% инициализация траектории движения UE1, м
dstep = 0.5;  % шаг сетки местоположений UE1 и UE2 по оси x и y, м
dstart1 = 0;  % начальная точка траектории движения UE1 по оси x, м
dstop1 = 30;  % конечная точка траектории движения UE1 по оси x, м
UE1_x_vec = dstart1:dstep:dstop1; % вектор точек траектории UE1 по оси x, м
UE1_x_vec_len=length(UE1_x_vec); % число точек траектории UE1 по оси x
UE1_z = 1.5;  % высота подвеса антенны UE1, м
UE1_y = 5;    % координата y, по которой двигается UE1, м
% инициализация траектории движения UE2, м
dstart2 = 30; % начальная точка траектории движения UE2, м
dstop2 = 0;   % конечная точка траектории движения UE1 по оси x, м
UE2_x_vec=dstart2:-dstep:dstop2; % вектор точек траектории UE2 по оси x, м
UE2_x_vec_len=length(UE2_x_vec); % число точек траектории UE2 по оси x
UE2_z = 1.5;  % высота подвеса антенны UE2, м
UE2_y_vec=5; % сценарий без разноса по оси y
% UE2_y = 0;    % координата y, по которой начинает двигаться UE2, м
% UE2_y_vec=UE2_y:dstep:UE2_y+10.0; % вектор точек траектории UE2 по оси y,м
UE2_y_vec_len=length(UE2_y_vec); % число точек траектории движения по оси y

%% формирование объектов диаграммообразования антенной решетки
% формирование объекта направляющего вектора для управления лучом АР
steervec_URA1 = phased.SteeringVector('SensorArray', URA1); % для URA1 gNB1
steervec_URA2 = phased.SteeringVector('SensorArray', URA2); % для URA2 gNB2
% формирование объекта вычисления коэффициента усиления антенной решетки
% при диаграммоообразовании в зависимости от весов направляющего вектора
gain_URA1 = phased.ArrayGain('SensorArray',      URA1, ...
                             'WeightsInputPort', true);
gain_URA2 = phased.ArrayGain('SensorArray',      URA2, ...
                             'WeightsInputPort', true);

% сетка азимутов и углов места для вычисления ДНА АР; для визуализации 
azimuth = -180:180;
elevation = -90:90;

% начальный азимут АР на gNB; для визуализации  
URA1_az = -270; % для URA1 на gNB1
URA2_az = 90;   % для URA2 на gNB2
% начальный угол места АР на gNB; нет наклона АР на gNB; для визуализации
URA1_elev = 0;  % для URA1 на gNB1
URA2_elev = 0;  % для URA2 на gNB2

ell_len = 0.1;  % начальное значение объекта эллипсоида местоположения UE
% масштабирующие коэффициенты для 3D объектов
STL_scale = .03;
AP_scale = 3;
use_image_scale = true; % предварительное масштабирование осей
if record_video == 1
    Filename = sprintf('lab_%s', datestr(now,'mm-dd-yyyy HH-MM'));
    recVideo = VideoWriter(Filename); % создание видеофайла для записи 
    recVideo.FrameRate = 10;          % частота кадров, достаточно 5 - 10
    open(recVideo);                   % открытие видеофайла для записи 
end

%% основной цикл имитационного моделирования
URA1UE1_SIR_vec = zeros(1,UE2_x_vec_len);   % вектор SIR для UE1
URA2UE2_SIR_vec = zeros(1,UE2_x_vec_len);   % вектор SIR для UE2
URA1UE1_SIR_mtx = zeros(UE2_y_vec_len, UE2_x_vec_len); % матрица SIR UE1
URA2UE2_SIR_mtx = zeros(UE2_y_vec_len, UE2_x_vec_len); % матрица SIR UE2
for UE2_y_idx=1:UE2_y_vec_len % цикл по координате y
    for d=1:UE2_x_vec_len % цикл по координате x
        if UE1_x_vec_len==1
            UE1_x = UE1_x_vec;
        else
            UE1_x = UE1_x_vec(d);
        end
        UE2_x = UE2_x_vec(d);
        UE2_y = UE2_y_vec(UE2_y_idx);
        % текущие координаты UE
        UE1_pos = [UE1_x; UE1_y; UE1_z];
        UE2_pos = [UE2_x; UE2_y; UE2_z];
        ell_len = 0.1;

        % оценка расстояния для вычисления потерь и углов между URA и UE
        [URA1_UE1_range, URA1_UE1_angle] = rangeangle(URA1_pos,UE1_pos); 
        [URA2_UE2_range, URA2_UE2_angle] = rangeangle(URA2_pos,UE2_pos); 
        [URA1_UE2_range, URA1_UE2_angle] = rangeangle(URA1_pos,UE2_pos); 
        [URA2_UE1_range, URA2_UE1_angle] = rangeangle(URA2_pos,UE1_pos); 

        % азитмуты и углы места от URA до UE для направляющего вектора
        URA1_UE1_az=URA1_UE1_angle(1)+180; URA1_UE1_el=-URA1_UE1_angle(2); 
        URA2_UE2_az=URA2_UE2_angle(1)+180; URA2_UE2_el=-URA2_UE2_angle(2); 
        URA1_UE2_az=URA1_UE2_angle(1)+180; URA1_UE2_el=-URA1_UE2_angle(2); 
        URA2_UE1_az=URA2_UE1_angle(1)+180; URA2_UE1_el=-URA2_UE1_angle(2); 
%         URA1_UE1_az=getAzimuth(URA1_pos,UE1_pos); URA1_UE1_el=-getElevation(URA1_pos,UE1_pos);
%         URA2_UE2_az=getAzimuth(URA2_pos,UE2_pos); URA2_UE2_el=-getElevation(URA2_pos,UE2_pos);
%         URA1_UE2_az=getAzimuth(URA1_pos,UE2_pos); URA1_UE2_el=-getElevation(URA1_pos,UE2_pos);
%         URA2_UE1_az=getAzimuth(URA2_pos,UE1_pos); URA2_UE1_el=-getElevation(URA2_pos,UE1_pos);

        % вычисление направляющих векторов SOI для АР gNB1 и gNB2
        sv_URA1 = steervec_URA1(fc,[URA1_UE1_az;URA1_UE1_el]);  % gNB1-UE1 
        sv_URA2 = steervec_URA2(fc,[URA2_UE2_az;URA2_UE2_el]);  % gNB2-UE2                 

        % вычисление КУ для оценки отношения сигнал/помеха SIR
        gain_URA1_UE1 = gain_URA1(fc,...
            [URA1_UE1_az; URA1_UE1_el], sv_URA1); % gNB1-UE1 SOI
        gain_URA2_UE2 = gain_URA2(fc,...
            [URA2_UE2_az; URA2_UE2_el], sv_URA2); % gNB2-UE2 SOI
        gain_URA1_UE2 = gain_URA1(fc,...
            [URA1_UE2_az; URA1_UE2_el], sv_URA1); % gNB1-UE2 SNOI
        gain_URA2_UE1 = gain_URA2(fc,...
            [URA2_UE1_az; URA2_UE1_el], sv_URA2); % gNB2-UE1 SNOI
        
        % заполнение вектора SIR для UE1
        URA1UE1_SIR_vec(d)=(gain_URA1_UE1 - fspl(URA1_UE1_range, lambda)) ...
                            - (gain_URA2_UE1 - fspl(URA2_UE1_range, lambda));
        % заполнение вектора SIR для UE2
        URA2UE2_SIR_vec(d)=(gain_URA2_UE2 - fspl(URA2_UE2_range, lambda))...
                            - (gain_URA1_UE2 - fspl(URA1_UE2_range, lambda));

        fprintf('Местоположение UE: y = %4.1f м, x = %4.1f м\n', UE2_y, UE2_x);
        fprintf('КУ в радиолинии SOI UE1-gNB1: = %5.1f дБи\n', gain_URA1_UE1);
        fprintf('КУ в радиолинии SNOI UE1-gNB2 = %5.1f дБи\n', gain_URA2_UE1);
        fprintf('КУ в радиолинии SOI  UE2-gNB2 = %5.1f дБи\n', gain_URA2_UE2);
        fprintf('КУ в радиолинии SNOI UE2-gNB1 = %5.1f дБи\n', gain_URA1_UE2);

        if plot_3D == 1
            % вычисление ДНА для двух URA по направляющему вектору
            [PAT_1, AZM_1, ELEV_1] = pattern(URA1, fc, azimuth, elevation, ...
                                            'CoordinateSystem', 'polar', ...
                                            'PropagationSpeed', c, ...
                                            'Type', 'powerdb', ...
                                            'Weights', sv_URA1);
            [PAT_2, AZM_2, ELEV_2] = pattern(URA2, fc, azimuth, elevation, ...
                                            'CoordinateSystem', 'polar', ...
                                            'PropagationSpeed', c, ...
                                            'Type', 'powerdb', ...
                                            'Weights', sv_URA2);

            % перевод дБ в линейные единицы для 3D визуализации
            PAT_1 = db2pow(PAT_1);
            PAT_2 = db2pow(PAT_2);

            % оценка полуплоскости yz (допустммы значения 1 и -1) командой
            % dot(q - p, n), где q - координата z URA, p - координата z UE,
            % 1 - нормаль
            plane_yz1 = (dot(URA1_pos(3) - UE1_pos(3), 1) >= 0) * 2 - 1;
            plane_yz2 = (dot(URA2_pos(3) - UE2_pos(3), 1) >= 0) * 2 - 1;

            % преобразование точек в сферическую систему координат (СК)
            [phi_1, theta_1] = meshgrid(AZM_1, ELEV_1);
            [phi_2, theta_2] = meshgrid(AZM_2, ELEV_2);
            rho_1 = PAT_1.*AP_scale;
            rho_2 = PAT_2.*AP_scale;
            x1 = rho_1.*sind(phi_1).*cosd(theta_1);
            y1 = rho_1.*sind(phi_1).*sind(theta_1);
            z1 = rho_1.*cosd(phi_1);
            x2 = rho_2.*sind(phi_2).*cosd(theta_2);
            y2 = rho_2.*sind(phi_2).*sind(theta_2);
            z2 = rho_2.*cosd(phi_2);

            % перемещение и поворот объектов в 3D
            x1_r = x1;
            y1_r = y1;
            z1_r = z1;
            [x1_r, y1_r, z1_r] = rotate_points(x1_r, y1_r, z1_r, 90, 'x');
            [x1_r, y1_r, z1_r] = rotate_points(x1_r, y1_r, z1_r, URA1_elev*plane_yz1, 'y');
            [x1_r, y1_r, z1_r] = rotate_points(x1_r, y1_r, z1_r, URA1_az, 'z');
            [x1_r, y1_r, z1_r] = move_points(x1_r, y1_r, z1_r, URA1_pos);

            x2_r = x2;
            y2_r = y2;
            z2_r = z2;
            [x2_r, y2_r, z2_r] = rotate_points(x2_r, y2_r, z2_r, 90, 'x');
            [x2_r, y2_r, z2_r] = rotate_points(x2_r, y2_r, z2_r, URA2_elev*plane_yz2, 'y');
            [x2_r, y2_r, z2_r] = rotate_points(x2_r, y2_r, z2_r, URA2_az, 'z');
            [x2_r, y2_r, z2_r] = move_points(x2_r, y2_r, z2_r, URA2_pos);

            % начало визуализации
            figh = figure;
            clf;

            % формирование эллипсоида местоположения UE
            [UE1_x_ell, UE1_y_ell, UE1_z_ell] = ellipsoid(UE1_pos(1), UE1_pos(2), UE1_pos(3), ...
                ell_len, ell_len, ell_len); % для UE1
            [UE2_x_ell, UE2_y_ell, UE2_ell] = ellipsoid(UE2_pos(1), UE2_pos(2), UE2_pos(3), ...
                ell_len, ell_len, ell_len); % для UE2 

            % отображение устройств UE и их соответствующих обозначений
            m1 = mesh(x1_r, y1_r, z1_r, rho_1);
            t1 = text(UE1_pos(1)+0.2, UE1_pos(2), 0, '\bfUE1', 'Rotation', 0);
            hold on;
            m2 = mesh(x2_r, y2_r, z2_r, rho_2);
            t2 = text(UE2_pos(1)+0.2, UE2_pos(2), 0, '\bfUE2', 'Rotation', 0);
            m3 = surf(UE1_x_ell, UE1_y_ell, UE1_z_ell);
            t3 = text(URA1_pos(1)+0.2, URA1_pos(2), 0, '\bfURA1', 'Rotation', 0);
            m4 = surf(UE2_x_ell, UE2_y_ell, UE2_ell);
            t4 = text(URA2_pos(1)+0.2, URA2_pos(2), 0, '\bfURA2', 'Rotation', 0);

            % установка границ графиков и обозначений осей
            maxX = max(max(x1, x2));
            maxY = max(max(y1, y2));
            maxZ = max(max(z1, z2));
            maxTotal = 2*max([maxX maxY maxZ]);
            xlim([-2*maxTotal 2*maxTotal]);
            ylim([-2*maxTotal 2*maxTotal]);
            zlim([-2*maxTotal 2*maxTotal]);
            xlabel("x"); ylabel("y"); zlabel("z");

            % отображение 3D объекта светофора и антенной решетки URA на нем 
            stl = stlread('trafficlight.stl');
            r = STL_scale; % масштабирующий коэффициент
            V = stl.vertices*r;
            [vx, vy, vz] = rotate_points(V(:,1), V(:,2), V(:,3), 90, 'x');
            [vx, vy, vz] = rotate_points(vx, vy, vz, 180, 'z');
            [vx1, vy1, vz1] = move_points(vx, vy, vz, URA1_pos.*STL_scale*1111);
            V1 = [vx1, vy1, vz1];
            p = patch(stl,'FaceColor',[0.8 0.8 1.0], ...
                'Vertices', V1*r, ...
                'EdgeColor', 'none', ...
                'FaceColor', [0 0 0], ...
                'FaceLighting', 'gouraud', ...
                'AmbientStrength', 0.15, ...
                'EdgeAlpha', 0);
            camlight('headlight');
            material('dull');

            st2 = stlread('trafficlight.stl');
            r = STL_scale; % scaling factor
            V = st2.vertices * r;
            [vx, vy, vz] = rotate_points(V(:,1), V(:,2), V(:,3), 90, 'x');
            [vx, vy, vz] = rotate_points(vx, vy, vz, 180, 'z');
            [vx2, vy2, vz2] = move_points(vx, vy, vz, URA2_pos.*STL_scale.*1111);
            V2 = [vx2, vy2, vz2];
            p = patch(st2,'FaceColor',[0.8 0.8 1.0], ...
                'Vertices', V2*r, ...
                'EdgeColor', 'none', ...
                'FaceColor', [0 0 0], ...
                'FaceLighting', 'gouraud', ...
                'AmbientStrength', 0.15, ...
                'EdgeAlpha', 0);

            % фиксация масштаба осей СК и установка угла наблюдения
            if use_image_scale
                axis('image');
            else
                xlim([dstart1*1.2, dstop1*1.2]); 
                ylim([dstart1*1.2, dstop1*1.2]);
                zlim([-1, 3]);
            end
            view([110 25]);

            % построение направляющей соединительной линии SOI между URA1 gNB1 и UE1
            plot3([URA1_pos(1), UE1_pos(1)], ...
                  [URA1_pos(2), UE1_pos(2)], ...
                  [URA1_pos(3), UE1_pos(3)], ...
                  'Color', 'g', 'LineWidth', 1.5);
            [~, a1] = rangeangle(UE1_pos, URA1_pos);
            [~, a2] = rangeangle(UE1_pos, URA1_pos);
            URA1UE1_angles = a2-a1;
            % построение направляющей соединительной линии SNOI между URA1 gNB1 и UE2
            plot3([URA1_pos(1), UE2_pos(1)], ...
                  [URA1_pos(2), UE2_pos(2)], ...
                  [URA1_pos(3), UE2_pos(3)], ...
                  '--', 'Color', 'r', 'LineWidth', 1.5);
            [~, a1] = rangeangle(UE1_pos, URA1_pos);
            [~, a2] = rangeangle(UE2_pos, URA1_pos);
            URA1UE2_angles = a2-a1;
            % построение направляющей соединительной линии SNOI между URA2 gNB2 и UE1
            plot3([URA2_pos(1), UE1_pos(1)], ...
                  [URA2_pos(2), UE1_pos(2)], ...
                  [URA2_pos(3), UE1_pos(3)], ...
                  '--', 'Color', 'r', 'LineWidth', 1.5);
            [~, a1] = rangeangle(UE2_pos, URA2_pos);
            [~, a2] = rangeangle(UE1_pos, URA2_pos);
            URA2UE1_angles = (a2-a1);
            % построение направляющей соединительной линии между SOI URA2 gNB2 и UE2
            plot3([URA2_pos(1), UE2_pos(1)], ...
                  [URA2_pos(2), UE2_pos(2)], ...
                  [URA2_pos(3), UE2_pos(3)], ...
                  'Color', 'g', 'LineWidth', 1.5);
            [~, a1] = rangeangle(UE2_pos, URA2_pos);
            [~, a2] = rangeangle(UE2_pos, URA2_pos);
            URA2UE2_angles = a2-a1;

            % проверка того, что мы не пересекли диапазон углов -180:180
            angle_list = [URA1UE1_angles, URA2UE1_angles, ...
                          URA1UE2_angles, URA2UE2_angles];
            for i=1:length(angle_list(1,:))
                angle_list(1, i) = wrapTo180(angle_list(1, i));
            end

            % отображение углов, соответствующих направлениям SOI и SNOI
            URA1UE1_angles = angle_list(:,1); % SOI  gNB1 URA1 - UE1
            URA2UE1_angles = angle_list(:,2); % SNOI gNB2 URA2 - UE1
            URA1UE2_angles = angle_list(:,3); % SNOI gNB1 URA1 - UE2
            URA2UE2_angles = angle_list(:,4); % SOI  gNB2 URA2 - UE2

            URA1UE1_angle_str = sprintf('%d^o', round(URA1UE1_angles(1)));
            tURA1UE1_angle = text(URA1_pos(1)+(UE1_pos(1)-URA1_pos(1))/2+0.1,...
                                  URA1_pos(2)+(UE1_pos(2)+URA1_pos(2))/2,...
                                  0,...
                                  URA1UE1_angle_str, 'Rotation', 0);
            URA1UE2_angle_str = sprintf('%d^o', round(URA1UE2_angles(1)));
            tURA1UE2_angle = text(URA1_pos(1)+(UE2_pos(1)-URA1_pos(1))/2+0.1,...
                                  URA1_pos(2)+(UE2_pos(2)-URA1_pos(2))/2,...
                                  0,...
                                  URA1UE2_angle_str, 'Rotation', 0);
            URA2UE1_angle_str = sprintf('%d^o', round(URA2UE1_angles(1)));
            tURA2UE1_angle = text(URA2_pos(1)+(UE1_pos(1)-URA2_pos(1))/2+0.2,...
                                  URA2_pos(2)+(UE1_pos(2)-URA2_pos(2))/2,...
                                  0,...
                                  URA2UE1_angle_str, 'Rotation', 0);
            URA2UE2_angle_str = sprintf('%d^o', round(URA2UE2_angles(1)));
            tURA2UE2_angle = text(URA2_pos(1)+(UE2_pos(1)-URA2_pos(1))/2+0.2, ...
                                  URA2_pos(2)+(UE2_pos(2)-URA2_pos(2))/2,...
                                  0,...
                                  URA2UE2_angle_str, 'Rotation', 0);
        end % if plot_3D == 1
        if record_video == 1
            pause(0.001)                 % пауза и захват кадра
            frame = getframe(gcf);       % получение кадра
            writeVideo(recVideo, frame); % запись кадра
        end
    end % for d=1:dist_vec_len
    URA1UE1_SIR_mtx(UE2_y_idx, :) = URA1UE1_SIR_vec; % матрица SIR для UE1
    URA2UE2_SIR_mtx(UE2_y_idx, :) = URA2UE2_SIR_vec; % матрица SIR для UE2
end % for UE2_y_idx=1:UE2_y_vec_len
if record_video == 1
    close(recVideo);    % закрытие видеофайла
end

%% визуализация результатов имитационного моделирования
if length(UE2_y_vec)==1
    figure(1000); 
    plot(UE1_x_vec, URA1UE1_SIR_mtx(1,:),'r'); hold on;
    plot(UE2_x_vec, URA2UE2_SIR_mtx(1,:),'b'); grid on;
    legend('SIR: SOI gNB_{1}->UE_{1}, SNOI: gNB_{2}->UE_{1}',...
           'SIR: SOI gNB_{2}->UE_{2}, SNOI: gNB_{1}->UE_{2}');
    xlabel('x, м'); ylabel('SIR = P_{SOI, дБ} - P_{SNOI, дБ}, дБ');
    title('Мгновенное отношение сигнал/помеха SIR'); hold off;
elseif length(UE2_y_vec)>1
    SIR_thr=20; % пороговое отношение сигнал/помеха, дБ

    URA1UE1_SIR_mtx=URA1UE1_SIR_mtx(:,end:-1:1);
    figure(1001); surf(UE2_x_vec, UE2_y_vec, URA1UE1_SIR_mtx);
    xlim([UE2_x_vec(end) UE2_x_vec(1)]); ylim([UE2_y_vec(1), UE2_y_vec(end)]);
    xlabel('x, м'); ylabel('y, м'); colormap default; shading interp; 
    c1_1 = colorbar; set(get(c1_1,'label'),'string','Мгновенное отношение SIR, дБ'); 
    view(2); set(0,'DefaultFigureWindowStyle','normal');
    title('Карта SIR в радиолинии SOI gNB_{1}->UE_{1}');

    figure(1002);   
    URA1UE1_SIR_mtx_area = double(URA1UE1_SIR_mtx > SIR_thr);
    surf(UE2_x_vec, UE2_y_vec, URA1UE1_SIR_mtx_area);
    xlim([UE2_x_vec(end) UE2_x_vec(1)]); ylim([UE2_y_vec(1), UE2_y_vec(end)]);
    xlabel('x, м'); ylabel('y, м'); shading interp; colormap winter; 
    colorbar; view(2); set(0,'DefaultFigureWindowStyle','normal');
    title(['Область с SIR > ',num2str(SIR_thr), ' дБ в радиолинии SOI gNB_{1}->UE_{1}']);

    figure(1003);    
    surf(UE2_x_vec, UE2_y_vec, URA2UE2_SIR_mtx);
    xlim([UE2_x_vec(end) UE2_x_vec(1)]); ylim([UE2_y_vec(1), UE2_y_vec(end)]);
    xlabel('x, м'); ylabel('y, м'); colormap default; shading interp; 
    c2_2 = colorbar; set(get(c2_2,'label'),'string','Мгновенное отношение SIR, дБ'); 
    view(2); set(0,'DefaultFigureWindowStyle','normal');
    title('Карта SIR в радиолинии SOI gNB_{2}->UE_{2}');

    figure(1004);   
    URA2UE2_SIR_mtx_area = double(URA2UE2_SIR_mtx > SIR_thr);
    surf(UE2_x_vec, UE2_y_vec, URA2UE2_SIR_mtx_area);
    xlim([UE2_x_vec(end) UE2_x_vec(1)]); ylim([UE2_y_vec(1), UE2_y_vec(end)]);
    xlabel('x, м'); ylabel('y, м'); shading interp; colormap winter; 
    colorbar; view(2); set(0,'DefaultFigureWindowStyle','normal');
    title(['Область с SIR > ',num2str(SIR_thr), ' дБ в радиолинии SOI gNB_{2}->UE_{2}']);
end


%% вспомогательные функции
function azimuth = getAzimuth(posvec1, posvec2)
azimuth = atan2d((posvec2(2) - posvec1(2)),(posvec2(1) - posvec1(1)));
end

function elevation = getElevation(posvec1, posvec2)
dZ = posvec2(3) - posvec1(3);
dX_sq = (posvec2(1) - posvec1(1))^2;
dY_sq = (posvec2(2) - posvec1(2))^2;
dZ_sq = (posvec2(3) - posvec1(3))^2;
elevation = 90 - acosd(abs(dZ) / sqrt(dX_sq + dY_sq + dZ_sq));
end

function [x_m, y_m, z_m] = move_points(x0, y0, z0, move_vec)
x = move_vec(1,:); y = move_vec(2,:); z = move_vec(3,:);
% перемещение вектора
x_m = x0 + x;      y_m = y0 + y;      z_m = z0 + z;
end
% function [x_m, y_m, z_m] = move_points(vec0, move_vec)
% x0 = vec0(1,:);       y0 = vec0(2,:);      z0 = vec0(3,:);
% x = move_vec(1,:);    y = move_vec(2,:);   z = move_vec(3,:);
% % move vector    
% x_m = x0 + x;         y_m = y0 + y;        z_m = z0 + z;
% end

function [x_r, y_r, z_r] = rotate_points(x, y, z, angle_deg, axis)
% поворот вектора 
switch axis
    case 'x'
        angle = angle_deg*pi/180;
        x_r = x;
        y_r = cos(angle).*y + -sin(angle).*z;
        z_r = sin(angle).*y + cos(angle).*z;
    case 'y'
        angle = angle_deg*pi/180;
        x_r = cos(angle).*x + sin(angle).*z;
        y_r = y;
        z_r = -sin(angle).*x + cos(angle).*z;
    case 'z'
        angle = angle_deg*pi/180;
        x_r = cos(angle).*x + -sin(angle).*y;
        y_r = sin(angle).*x + cos(angle).*y;
        z_r = z;
    case 'xyz'
        angle = angle_deg.*pi./180; 
        roll = angle(1);
        pitch = angle(2);
        yaw = angle(3);
        x_r = x;
        y_r = y;
        z_r = z;
        
        y_r = cos(roll).*y_r + -sin(roll).*z_r;
        z_r = sin(roll).*y_r + cos(roll).*z_r;
        
        x_r = cos(pitch).*x_r + sin(pitch).*z_r;
        z_r = -sin(pitch).*x_r + cos(pitch).*z_r;
        
        x_r = cos(yaw).*x_r + -sin(yaw).*y_r;
        y_r = sin(yaw).*x_r + cos(yaw).*y_r;
    case 'zyx'
        angle = angle_deg.*pi./180; 
        roll = angle(1);
        pitch = angle(2);
        yaw = angle(3);
        x_r = x;
        y_r = y;
        z_r = z;
        x_r = cos(yaw).*x_r + -sin(yaw).*y_r;
        y_r = sin(yaw).*x_r + cos(yaw).*y_r;
        
        x_r = cos(pitch).*x_r + sin(pitch).*z_r;
        z_r = -sin(pitch).*x_r + cos(pitch).*z_r;
        
        y_r = cos(roll).*y_r + -sin(roll).*z_r;
        z_r = sin(roll).*y_r + cos(roll).*z_r;
    otherwise
        warning('Unexpected axis type')
    end
end