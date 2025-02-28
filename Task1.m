% Чтение данных из файла
data = load('skirk001.txt');

% Извлечение данных из столбцов
x = data(:, 1);
y1 = data(:, 2);
y2 = data(:, 3);

% Построение графиков
figure;

% Первый график
% Первый график с толстой линией
plot(x, y1, 'g-', 'LineWidth', 2, 'DisplayName', 'График x(t)');
hold on;

% Второй график с толстой линией
plot(x, y2, 'b--', 'LineWidth', 2, 'DisplayName', 'График v(t)');


% Настройка графика
xlabel('t');
xlim([0, 15]);
ylabel(' ');

title('Графики разности решений с шагом сетки h=0,01', 'Fontsize', 22);
legend show;
set(h_leg, 'FontSize', 34);
grid on;

hold off;

