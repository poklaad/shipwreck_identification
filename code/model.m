% Решается дифур 2-ого порядка: theta'' + nu*theta' + rec(theta,time) = waves(time, A) 
% A*cos(omega*time) - регулярное волнение




% ПАРАМЕТРЫ ДЛЯ ПОСТРОЕНИЯ ДИФФЕРЕНЦИАЛЬНОГО УРАВНЕНИЯ %%%%%%%%%%%%%%%%%%%

% y0_1, y0_2 - начальные условия для решения дифура
y0_1 = 0;
y0_2 = 0;

% nu - параметр из дифура при theta'
% A - высота волнения
% T - период волны, с. Для регулярного волнения
% omega - параметр из дифура в cos - частота волнения. Для регулярного волнения
nu = 0.01;
A = 0.05;
%T = 2*pi/0.64;
%omega = 2*pi/T;





% МОДЕЛИРОВАНИЕ РАЗВИТИЯ АВАРИИ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Колебания корабля в процессе затопления моделируются в цикле.

% alltime - длина отрезка моделирования
% start_time - параметр, определяющий момент начала моделирования.
% start_wreck - время начала развития аварии.
alltime = 4000;
start_time = 1;
start_wreck = 200;

% step - шаг времени, который определяет, как часто берутся точки для
% построения графиков. С его помощью можно управлять частотой сбора данных.
% По умолчанию равен 1 секунде.
step = 1;

% pose1 - начальный тип статической остойчивости коробля - всегда 1-ый непорежденный
% pose2 - конечный тип статической остойчивости коробля
pose1 = 1;
pose2 = 5;





% ФУНКЦИЯ ВОССТАНАВЛИВАЮЩЕГО МОМЕНТА %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a0, a1, a3, a5 - коэффициенты полинома пятого порядка, задающего диаграмму статической остойчивости коробля.
% Порядок типов остойчивости: 
% 	1-ый неповрежденный     1-ый		2-ой		3-ий		4-ый		5-ый

a0=[		0               0           -0.2		0       	-0.2 		0.07]; 
a1=[		0.64            0.25		0.64		-0.64		-0.64		-0.64]; 
a3=[		-0.1            -0.1		-0.1 		2.5 		2.5 		2.5]; 
a5=[		-0.07           -0.05		-0.07		-1.3		-1.3		-1.3];





% ВОЛНЕНИЕ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Файлы с волнением. Содержат углы волнового склона в градусах с шагом 1 секунда
file_id = fopen('ANG4.DAT','r'); % волновой склон 4 балла
%file_id = fopen('ANG5.DAT','r'); % волновой склон 5 баллов
%file_id = fopen('ANG6.DAT','r'); % волновой склон 6 баллов
%file_id = fopen('ANG7.DAT','r'); % волновой склон 7 баллов
%file_id = fopen('ANG8.DAT','r'); % волновой склон 8 баллов
%file_id = fopen('anglM.DAT','r'); % ветровое волнение
%file_id = fopen('anglS.DAT','r');% зыбь
%file_id = fopen('anglWW.DAT','r'); % смешанное волнение

% angle_waves - углы волнового склона в градусах с шагом 1 секунда
angle_waves = fscanf(file_id,'%f');

fclose(file_id);

% wv - величина волнения по оси ординат
% x_waves - моменты измерения волнения во времени
wv = [0];
%wv = waves (0, alltime, 1, 50);
x_waves = start_time:step:alltime;

% Вычисление высоты волны в каждый момент времени
for i = 2:length(angle_waves)
    wv(i) = tand(angle_waves(i-1));
end

% Масштабирование волнения. А - максимальная по модулю волна
[tmp, ind1] = max(abs(wv));
scale_wave = 0.05 / tmp;
wv = scale_wave * wv;

% Просмотр всего доступного волнения
x_waves = start_time : 1 : length(wv);
figure
plot(x_waves(1:9000), wv(1:9000), "k"); 
axis([2 9000 -0.05 0.05]) 
% ограничение промежутка волнения до требуемой длины
wv = wv(1:length(x_waves));


% Интерполирование волнения
fun_waves = @(new_x) interp1(x_waves, wv, new_x, 'spline');





% ИНИЦИАЛИЗАЦИЯ ПЕРЕМЕННЫХ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f1 - окно для фазового портрета колебаний корабля
% f2 - окно для колебаний за последние aver_inter секунд
% f3 - окно для отображения всех колебаний
% f4 - окно для квазистационарных колебаний
% f5 - окно для функции автокорреляции_1
% f6 - окно для функции автокорреляции_2
% f7 - окно для концентрических окружностей (функция concentric_circles)
%f1=figure;
%f2=figure;
f3=figure;
f4=figure;
%f5=figure;
%f6=figure;
%f7=figure;


% theta_data и d_theta_data - массивы для хранения функции колебаний корабля и ее производной соответственно
% time_data - массив для хранения моментов времени, в которых вычислена функция колебаний
theta_data = [];
d_theta_data = [];
time_data = [];

% expectation и dispersion - массивы для хранения выборочного среднего и дисперсии колебаний
expectation = [];
dispersion = [];


% circles_count - количество когнитивных спиралей
%circles_count = 30;


% r_count - количество корней полинома, задаюшего функцию восстанавливающего момента.
% old_wreck_type, new_wreck_type - вычисленные типы затопления на предыдущем и нынешнем шагах
r_count = 0;
old_wreck_type = 1;
new_wreck_type = 1;

% right_inflection_point, old_right_inflection_point - новое и старое
% значение функции восстанавливающего момента в правой точке перегиба.
% Используется для идентификации 1-ого типа затопления
old_right_inflection_point = 0;
right_inflection_point = 0;

% inflection_points_count - количетсво точек перегиба.
% changed_wreck_type - флаг для определения, изменился ли типа затопления.
% where_changed - точки, где изменился тип затопления
% wreck_types_through_modeling - массив с типами затопления, в которых находилось судно
inflection_points_count = 0;
changed_wreck_type = false;
where_changed = [1];
wreck_types_through_modeling = [2];

% функция автокорреляции и задержка для ее вычисления
autocor_theta = [];
autocor_d_theta = [];
delay = 200;




% ОПРЕДЕЛЕНИЕ КВАЗИСТАЦИОНАРНОСТИ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aver_inter - длина интервала, по которому производится осреднение.
% Примерно равен длительности 18-ти колебаний судна в секундах.
% quas_len - минимальная длина участка, на котором среднее и дисперсия
% должны быть стабильны, чтобы считать его квазистационарным
% is_quas_stat - проверка, находится ли сейчас судно в квазистационарном состоянии
aver_inter = 144;
quas_len = 80;
is_quas_stat = false;

% quas_starts - начала квазистационарных участков
% quas_ends - концы квазистационарных участков
quas_starts = [0];
quas_ends = [0];

% accept_mean - допустимое изменение выборочного среднего на квазистационарном участке
% accept_disp - допустимое изменение дисперсии на квазистационарном участке
accept_mean = 0.5;
accept_disp = 0.026;

% exp_direction - вычисленное направление среднего на квазистационарном участке
% new_exp_direction - направление среднего по последней точке на квазистационарном участке
% accept_dir - допустимое изменение направления среднего на квазистационарном участке
% len_dir_calc - длина участка для вычисления направления среднего на
% квазистационарном участке
direction = 0;
new_exp_direction = 0;
accept_dir = 0.4;
len_dir_calc = 30;





% ДЕТРЕНД %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% detrend_step - длина участка, на котором производится детренд
% detrend_theta_data - массив колебаний с детрендом
detrend_step = 24*step;
detrend_theta_data = [];





% МОДЕЛИРОВАНИЕ ПРОЦЕССА РАЗВИТИЯ АВАРИИ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time = start_time:step:alltime
    
    % ВЫЧИСЛЕНИЕ ПОЛОЖЕНИЯ СУДНА %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % proc - процент прошедшего времени. Используется для вычисления нынешнего значения коэффициентов диаграммы остойчивости
    if time <= start_wreck
        proc = 0;
    else
        proc = (time - start_wreck) / alltime;
    end
    %proc = 0;
    
    % recovery - нынешние коэффициенты диаграммы остойчивости.
    recovery = [a0(pose1) + ((a0(pose2) - a0(pose1)) * proc) 
                a1(pose1) + ((a1(pose2)-a1(pose1)) * proc) 
                0 
                a3(pose1) + ((a3(pose2)-a3(pose1)) * proc) 
                0 
                a5(pose1) + ((a5(pose2)-a5(pose1)) * proc) 
                ];
       
            
    % DoDt - вспомогательная функция (расположена в конце программы) для
    % вычисления дифура
    [time_dif,theta_dif] = ode45(@(t,y) DoDt(t, y, fun_waves, nu, recovery),[(time-step) time],[y0_1 y0_2]);

    
    
    
    
    % ЗАПИСЬ ПОЛОЖЕНИЯ КОРОБЛЯ В МОМЕНТ TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Каждому моменту времени соответствует одна точка на графике.
    % Точность/гладкость графика регулируется параметром step.
    % Чем меньше step, тем более гладкий график, однако уменьшение step
    % отрицательно влияет на время работы программы.
    theta_data(end + 1) = theta_dif(end, 1);
    d_theta_data(end + 1) = theta_dif(end, 2);
    time_data(end + 1) = time_dif(end);
    
    % При решении дифура на следующем шаге цикла будет использоваться его
    % состояние в конце нынешнего шага цикла
    y0_1 = theta_dif(end, 1);
    y0_2 = theta_dif(end, 2);
    
   
    
    
    
    % ОПРЕДЕЛЕНИЕ ТИПА ЗАТОПЛЕНИЯ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_wreck_type = new_wreck_type;
    if time == 1
        dif_recovery = polyder([a5(pose1) 0 a3(pose1) 0 a1(pose1) a0(pose1)]);
        dif_roots = roots(dif_recovery);
        old_right_inflection_point = polyval([a5(pose1) 0 a3(pose1) 0 a1(pose1) a0(pose1)], dif_roots(end));
    end
    
    % Вычисление особенностей (корней и точек перегиба) функции восстанавливающего момента
    
    % Корни
    recovery = flipud(recovery);
    r = roots(recovery);
    not_im = imag(r)==0;
    r = r(not_im);
    r_count = length(r);
    r = sort (r);
    
    % Точки перегиба
    dif_recovery = polyder(recovery);
    dif_roots = roots(dif_recovery);
    not_im_dif_roots = imag(dif_roots)==0;
    dif_roots = dif_roots(not_im_dif_roots);
    inflection_points_count = length(dif_roots);
    right_inflection_point = polyval(recovery, dif_roots(end));
    
    
    % Проверка на 5-ый тип затопления (5 корней, несимметричны)
    if r_count == 5 && abs(1 - abs(r(2))/abs(r(4))) > 0.1
        new_wreck_type = 6;
    else
        % Проверка на 4-ый тип затопления (3 корня, есть перегиб)
        if r_count == 3 && inflection_points_count == 4
            new_wreck_type = 5;
        else
            % Проверка на 3-ий тип затопления (5 корней, симметричны)
            if r_count == 5 && abs(1 - abs(r(2))/abs(r(4))) <= 0.1
                new_wreck_type = 4;
            else
                % Проверка на 2-ой тип затопления (3 корня, второй корень не в нуле)
                if r_count == 3 && r(2) ~= 0
                    new_wreck_type = 3;
                else
                    % Проверка на 1-ый тип затопления (правая точка перегиба опустилась)
                    if right_inflection_point < old_right_inflection_point
                        new_wreck_type = 2;
                    end
                end
            end
        end
    end       
    
    % Если изменился тип затопления - добавить точку
    if old_wreck_type ~= new_wreck_type
        changed_wreck_type = true;
        where_changed = [where_changed time];
        wreck_types_through_modeling = [wreck_types_through_modeling new_wreck_type];
    end
    
    
    
    

    % ВЫВОД ГРАФИКОВ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Определение количества точек для отображения на фазовом портрете (по
    % умолчанию равно длине интервала для осреднения)
    points_to_show = points_for_aver;
    
    % Выделение массивов для отрисовки
    theta_data_show = theta_data(points_to_show:end);
    d_theta_data_show = d_theta_data(points_to_show:end);
    time_data_show = time_data(points_to_show:end);
    expectation_show = expectation(points_to_show:end);
    dispersion_show = dispersion(points_to_show:end);
    % Поиск точки времени, до которой нужно отобразить график волнения
    [tmp, time_wv] = min(abs(x_waves - time_data(end)));
    x_waves_show = x_waves(points_to_show:time_wv);
    wv_show = wv(points_to_show:time_wv);
    
    % Вывод графика фазового портрета
    figure(f1);
    clf(f1);
    plot (theta_data_show, d_theta_data_show, "b")
    axis([-1 1 -1 1])
    drawnow 

    % Вывод графиков колебаний корабля, волнения, выборочного среднего и дисперсии
    figure(f2);
    clf(f2);
    hold on
    plot(time_data_show, theta_data_show, "b")
    plot(time_data_show, expectation_show, "r")
    plot(time_data_show, dispersion_show, "g")
    
    plot(x_waves_show, wv_show, "k");
    % Отображение точек, в которых изменился тип затопления
    changes_show = where_changed(where_changed >= time - aver_inter);
    if changed_wreck_type && ~isempty(changes_show)
        xline(changes_show);
    end
    x1 = max(time - aver_inter - 10, 0);
    x2 = time + 10;
    y1 = min(min(theta_data), -0.1);
    y2 = max(max(theta_data), 0.1);
    axis([x1 x2 y1 y2])
    hold off
end





% ДЕТРЕНД %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time = start_time:detrend_step:alltime - detrend_step
    detrend_section = detrend(theta_data(time:time + detrend_step - 1));
    detrend_theta_data = [detrend_theta_data detrend_section];
end
detrend_section = detrend(theta_data(time + detrend_step:end));
detrend_theta_data = [detrend_theta_data detrend_section];





% ВЫЧИСЛЕНИЕ СТАТИСТИЕСКИХ ХАРАКТЕРИСТИК %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time = start_time:step:alltime
    % ВЫЧИСЛЕНИЕ СРЕДНЕГО И ДИСПЕРСИИ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Осредняется последняя часть графика длины aver_inter секунд.
    
    % Поиск начальной точки для осреднения
    [tmp, points_for_aver] = min(abs(time - aver_inter - time_data));
    % Если еще не прошлого достаточно времени, осредняется весь график.
    points_for_aver = max(1, points_for_aver);
    expectation(end + 1) = mean(detrend_theta_data(points_for_aver:time / step));
    dispersion(end + 1) = var(detrend_theta_data(points_for_aver:time / step));
    
    
    
    
    
    % ПОИСК ИНТЕРВАЛОВ КВАЗИСТАЦИОНАРНОСТИ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %if time >= start_wreck
        % Если квазистационарный участок закончился, начинаем его
        if ~is_quas_stat
            is_quas_stat = true;
            quas_starts = [quas_starts time];
        end
        % Направление среднего обновляется в первые len_dir_calc секунд 
        % квазистационарного участка. После их прошествия оно фикисруется
        if time - quas_starts(end) * step <= len_dir_calc
            mean_exp =  mean(expectation(quas_starts(end):end));
            ind = quas_starts(end)*step + (length(expectation) - quas_starts(end)) * step / 2;
            direction = (2*mean_exp - expectation(quas_starts(end)))/(time - quas_starts(end) * step + 1);
            direction = 0;
            
            % Поиск экстремумов колебаний
            maxs = islocalmax(expectation);  
            mins = islocalmin(expectation);    
            % Индексы первых максимумов и минимумов в массиве графика колебаний 
            if length(maxs) ~= 1
                quas_maxs_ind = find(maxs);
            else 
                quas_maxs_ind = 1;
            end
            if length(mins) ~= 1
                quas_mins_ind = find(mins);
            else 
                quas_mins_ind = 1;
            end
            max_ind = find(quas_maxs_ind > quas_starts(end));
            quas_maxs_ind = quas_maxs_ind(max_ind);
            min_ind = find(quas_mins_ind > quas_starts(end));
            quas_mins_ind = quas_mins_ind(min_ind);
            mean_max = mean(expectation(quas_maxs_ind));
            mean_min = mean(expectation(quas_mins_ind));
            
            add_comp_max = mean_max - ind * direction;
            add_comp_min = mean_min - ind * direction;
        end
        
         % Поиск экстремумов колебаний
        maxs = islocalmax(expectation);  
        mins = islocalmin(expectation);    
        % Индексы первых максимумов и минимумов в массиве графика колебаний 
        last_maxs_ind = find(maxs,5,'last');
        last_mins_ind = find(mins,5,'last');
        
        maxs_count = 0;
        mins_count = 0;
        if length(expectation) > 1
            for i=1:length(last_maxs_ind)
                if abs(expectation(last_maxs_ind(i))) > abs((direction * last_maxs_ind(i) + add_comp_max) * (1 + accept_dir)) || abs(expectation(last_maxs_ind(i))) < abs((direction * last_maxs_ind(i) + add_comp_max) * (1 - accept_dir))
                    maxs_count = maxs_count + 1;
                end
            end
            for i=1:length(last_mins_ind)
                if abs(expectation(last_mins_ind(i))) > abs((direction * last_mins_ind(i) + add_comp_min) * (1 + accept_dir)) || abs(expectation(last_mins_ind(i))) < abs((direction * last_mins_ind(i) + add_comp_min) * (1 - accept_dir))
                    mins_count = mins_count + 1;
                end
            end
        end
        
        % Если направление среднего изменилось больше, чем в accept_dir
        % раз, участок неквазистационарный
        if (maxs_count == 5 || mins_count == 5) && time - quas_starts(end) * step >= len_dir_calc
            % Участок неквазистационарный, заканчиваем его.
            is_quas_stat = false;
            % Удаляются последние 5 колебаний
            quas_end_time = min(last_mins_ind(1), last_maxs_ind(1));
            % Участок добавляется, если он больше минимальной длины
            if quas_end_time - quas_starts(end) >= quas_len
                % Длинные квазистационарные участки делятся на участки
                % длинны 100 секунд
                for quas_parts = 100:100/step:quas_end_time - quas_starts(end)
                    quas_ends = [quas_ends (quas_starts(end) + 100/step)];
                    quas_starts = [quas_starts (quas_starts(end) + 100/step)];
                end
                quas_starts(end) = [];
                quas_ends(end) = quas_end_time;
            else 
                quas_starts = quas_starts(1:length(quas_starts)-1);
            end
            % Начинаем новый участок
        end
    %end
end

% Если до конца длился квазистационарный участок, завершаем его
if is_quas_stat 
    is_quas_stat = false;
    if time - step - quas_starts(end) >= quas_len
        quas_ends = [quas_ends (time - step)];
    else 
        quas_starts = quas_starts(1:length(quas_starts)-1);
    end
end
disp('End of staticstics calculating')





% ВЫВОД КВАЗИСТАЦИОНАРНЫХ УЧАСТКОВ
figure(f4);
% piece_len - общая длина квазистационарных участков
piece_len = 0;
quas_theta = [];
quas_d_theta = [];
quas_disp = [];
for i = 2:length(quas_starts)
    piece_len = piece_len + quas_ends(i) - quas_starts(i) + 1;
    quas_theta = [quas_theta detrend_theta_data(quas_starts(i):quas_ends(i))];
    quas_d_theta = [quas_d_theta d_theta_data(quas_starts(i):quas_ends(i))];
    quas_disp = [quas_disp dispersion(quas_starts(i):quas_ends(i))];
    
end
quas_time = (1:step:piece_len);
plot (quas_time, quas_theta);
    
    
    
    
    
% ВЫЧИСЛЕНИЕ АВТОКОРРЕЛЯЦИИ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Автокорреляция вычисляется на квазистацонарных участках на каждом 
% шаге моделирования с задержкой delay
for i = 2:length(quas_starts)
    detrend_quas = detrend(theta_data(quas_starts(i):quas_ends(i)),1);
    [cor_now, del_lags] = xcorr(detrend_quas, delay, 'normalized');
    %[d_cor_now, del_lags] = xcorr(d_theta_data(quas_starts(i):quas_ends(i)), delay);
    cor_now = cor_now(round(length(cor_now)/2):end);
    peaks = islocalmax(cor_now);  
    % Индексы первых максимумов и минимумов в массиве графика колебаний 
    peaks_ind = find(peaks);
    peaks = cor_now(peaks);
    cut = peaks_ind(end);
    for j = 1:length(peaks)-1
        if peaks(j) < peaks(j + 1)
            cut = peaks_ind(j);
            cut = delay;
            break;
        end
    end
    cor_now = cor_now(1:cut);
    cor_now = [cor_now zeros(1, delay - length(cor_now))];
    autocor_theta = [autocor_theta, cor_now.'];
end
    




% ВЫВОД ПОЛНОГО ГРАФИКА КОЛЕБАНИЙ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(f3);
hold on
plot(time_data , theta_data, "b")
%plot(time_data , detrend_theta_data, "b")
plot(time_data , expectation, "r")
plot(time_data , dispersion, "g")
% Выделение цветом квазистационарных участков
for i = 2:length(quas_starts)
    plot(time_data(quas_starts(i):quas_ends(i)) , theta_data(quas_starts(i):quas_ends(i)))
    %plot(time_data(quas_starts(i):quas_ends(i)) , detrend_theta_data(quas_starts(i):quas_ends(i)))
end
plot(x_waves , wv, "k");
% Отображение точки, в которой изменился тип затопления
if changed_wreck_type
    xline(where_changed);
end
hold off





% ВЫВОД АВТОКОРРЕЛЯЦИИ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Автокорреляция для каждого квазистационарного участка с задержкой delay

for j = 2:length(quas_starts)
    figure;
    hold on
    % отрезание нулевой части
    show_ind = find (autocor_theta(:,j-1) ~= 0,1,'last');
    autocor_show = autocor_theta(1:show_ind,j-1);
    del_lags = (0:length(autocor_show)-1);
    plot(del_lags, autocor_show, "b")
    hold off
end





% ЗАПИСЬ В ФАЙЛ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_autocor = fopen('autocor.txt','a');
file_wreck_type = fopen('wreck_type.txt','a');
% Минимальная длина функции автокорреляции
minimal_autocor_length = 30;
for i = 2:length(quas_starts)
    
    % Выделение функции автокорреляции заданной длины
    show_ind = find (autocor_theta(:,i-1) ~= 0,1,'last');
    autocor_show = autocor_theta(1:show_ind,i-1);
    del_lags = (0:length(autocor_show)-1);
    
    if length(autocor_show) > minimal_autocor_length
        autocor_show = autocor_show(1:minimal_autocor_length);

        
        % Поиск типов затопления на данном квазистационарном участке
        ind_changes_before_end = find(where_changed <= quas_ends(i));
        ind_changes_after_start = find(where_changed >= quas_starts(i) & where_changed <= quas_ends(i));
        ind_changes_before_end = ind_changes_before_end(length(ind_changes_before_end) - length(ind_changes_after_start):end);
        wrecks_write = wreck_types_through_modeling(ind_changes_before_end)-1;
        % Если данному участку соответствует несколько типов
        % затопления, каждый из них записывается отдельно, и
        % для каждого из них отдельно записывается функция
        % автокорреляции
        for type_count = 1:length(wrecks_write)
            fprintf(file_autocor,'%.6f ',autocor_show);
            fprintf(file_autocor,'\n');
            fprintf(file_wreck_type,'%i ',wrecks_write(type_count));
            fprintf(file_wreck_type,'\n');
        end
    end
end
fclose(file_autocor);
fclose(file_wreck_type);






% ДИФФЕРЕНЦИАЛЬНОЕ УРАВНЕНИЕ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DoDt реализует дифур 2-ого порядка в виде системы с заменой переменных:
% {d th(1)/ dt = th(2)
% {d th(2)/ dt = A*cos(omega*time) - nu*th(2) - recovery(th(1))

function DthetaDtime = DoDt(t, th, fun_waves, nu, recovery)
    DthetaDtime = [th(2);fun_waves(t) - nu*th(2) - (recovery(1)+recovery(2)*th(1)+recovery(4)*th(1).^3+recovery(6)*th(1).^5)];
end