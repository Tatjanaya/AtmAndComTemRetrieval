% fileNameInput TIGR���������ļ�
% fileNameOutput �����tape5�ļ�
% fltPath Sentinel-3 SLSTR ������Ӧ����
% kind ����Ĵ�������
% angle �Ƕ�
function tape5Result = createTape5(fileNameInput, fileNameOutput, fltPath, kind, angle)
    % ��946�����������ļ� 1-236-�ȴ� 237-375-��γ���ļ� 376-494-��γ�ȶ��� 495-573-�����ļ�
    % 574-946-���ض���
    fidin = fopen(fileNameInput, 'r');
    fidout = fopen(fileNameOutput, 'a');
    for i = 1:946
       % �жϼ��ڣ�Ĭ���ȴ�
       season = 1;
       if i <= 236
           season = 1;
       elseif i <= 375 && i >= 237
           season = 2;
       elseif i <= 494 && i >= 376
           season = 3;
       elseif i <= 573 && i >= 495
           season = 4;
       else
           season = 5;
       end
       
       % �����ַ���
       seasonStr = strcat(season, '    ');
       seasonStrLink = sprintf('%s%s%s%s', seasonStr, seasonStr, seasonStr, seasonStr);
       lines = textscan(fidin,'%f%f%f%f%f',40,'Headerlines',(i - 1) * 43 + 2);
       
       % AAH Cxxx
       AAHseason = season + 1;
       AAHseasonStr = sprintf('%s%d%d%d', 'AAH C', AAHseason, AAHseason, AAHseason);
       
       % card5 IRPT
       IRPT = 1;
       if i == 946
           IRPT = 0;
       end
       
       % ��ȡcell�и��в�ת��Ϊmat
       altitude = cell2mat(lines(1));
       pressure = cell2mat(lines(2));
       temperature = cell2mat(lines(3));
       waterVapor = cell2mat(lines(4));
       ozone = cell2mat(lines(5));
       
       % row col
       [row, col] = size(altitude);
       
       % �жϴ�ʱ���������ִ������� kind = 1��͸���� kind = 2������ kind = 3������
       if kind == 1
           fprintf(fid, sprintf('%s%s%s', 'TMF 7    2    0   -1    0    0    ', seasonStrLink, '0    1    0   0.001    0.0   !card1'));
           fprintf(fid, '\r\n');
           fprintf(fid, 'FFF  8 0.0   365.000                    01 F T T       0.000      0.00     0.000     0.000     0.000         0   !card1a');
           fprintf(fid, '\r\n');
           fprintf(fid, fltPath);
           fprintf(fid, '\r\n');
           fprintf(fid, '    1    0    0    0    0    0    23.000     0.000     0.000     0.000     0.500   !card2');
           fprintf(fid, '\r\n');
           fprintf(fid, '   40    0    0                           0.0    0     1.000    28.964  !card2c');
           fprintf(fid, '\r\n');
           % ˳�������Ǹ߳� ѹǿ ���� ˮ������ ������̼��0������ +AAH C222  2��ʾ����������γ���ļ�
           for j = 1:row
               fprintf(fid, sprintf('%s%4.3f%s%s%s%s%s%s%s%s%s%s%s', '   ', altitude(row, 1), ''));
           end
       elseif kind == 2
           fprintf(fid, sprintf('%s%s%s', 'TMF 7    2    1   -1    0    0    ', seasonStrLink, '0    1    0   0.001    0.0   !card1'));
           fprintf(fid, '\r\n');
           fprintf(fid, 'FFF  8 0.0   365.000                    01 F T T       0.000      0.00     0.000     0.000     0.000         0   !card1a');
           fprintf(fid, '\r\n');
           fprintf(fid, fltPath);
           fprintf(fid, '\r\n');
           fprintf(fid, '    1    0    0    0    0    0    23.000     0.000     0.000     0.000     0.500   !card2');
           fprintf(fid, '\r\n');
           fprintf(fid, '   40    0    0                           0.0    0     1.000    28.964  !card2c');
           fprintf(fid, '\r\n');
           % ˳�������Ǹ߳� ѹǿ ���� ˮ������ ������̼��0������ +AAH C222  2��ʾ����������γ���ļ�
           for j = 1:row
               
           end
       else
           fprintf(fid, sprintf('%s%s%s', 'TMF 7    2    1   -1    0    0    ', seasonStrLink, '0    1    0   0.001    0.0   !card1'));
           fprintf(fid, '\r\n');
           fprintf(fid, 'FFF  8 0.0   365.000                    01 F T T       0.000      0.00     0.000     0.000     0.000         0   !card1a');
           fprintf(fid, '\r\n');
           fprintf(fid, fltPath);
           fprintf(fid, '\r\n');
           fprintf(fid, '    1    0    0    0    0    0    23.000     0.000     0.000     0.000     0.500   !card2');
           fprintf(fid, '\r\n');
           fprintf(fid, '   40    0    0                           0.0    0     1.000    28.964  !card2c');
           fprintf(fid, '\r\n');
           % ˳�������Ǹ߳� ѹǿ ���� ˮ������ ������̼��0������ +AAH C222  2��ʾ����������γ���ļ�
           for j = 1:row
               
           end
       end
    end
end