classdef viewer < handle
% viewer 三维轨道绘图对象：按场景初始化背景，再用 touchSat 添加卫星轨道
%
% 场景（构造参数）：
%   'Earth'    场景一：3D 地球（J2000，默认）
%   'Moon'     场景二：3D 月球（月心惯性系）
%   'Synodic'  场景三：地月会合坐标系（旋转系；归一化单位，长度 = 地月距离；
%              原点为地月质心，地球 [-mu,0,0]，月球 [1-mu,0,0] 静止不动；
%              地球、月球与 L1~L5 均以等大小散点绘制（淡蓝/灰/深灰））
%
% 用法：
%   v = viewer('Earth');                              % 创建 3D 地球图
%   v = viewer('Moon');                               % 创建 3D 月球图
%   v = viewer('Synodic');                            % 创建地月会合系图
%   v.touchSat('../data/output/ISS_trajectory.txt');  % 添加一条卫星轨道
%   v.touchSat('../data/output/preciseEph_C06.txt');  % 再添加一条（可多次调用）
%   v.touchSat(fn, 'C06');                            % 可指定卫星名字，缺省为 sat1, sat2, ...
%   v.killSat('C06');                                 % 从图中删除卫星
%   v.zoomTo('L2');                                   % 放大观察局部：中心移到 L2（Synodic 场景）
%   v.zoomTo('L2', 0.05);                             % 自定义窗口半宽（归一化单位）
%   v.zoomTo([x, y, z], h);                           % 中心移到任意点（Synodic 为归一化坐标）
%   v.zoomOut();                                      % 恢复全景
%
% 每添加/删除一条卫星，图例自动重建，每条卫星显示一行：卫星名
%
% 星历文件格式：
%   8 列：MJD日, MJD秒, X, Y, Z, VX, VY, VZ（km, km/s）
%   或 7 列：MJD, X, Y, Z, VX, VY, VZ
% 坐标单位均为 km；'Synodic' 场景输入为地月质心为原点的会合系坐标（km），
% 内部自动除以地月距离归一化后绘制

    properties
        scene = 'Earth' % 场景：'Earth' / 'Moon' / 'Synodic'
        fig             % 图形窗口
        ax              % 坐标轴
        earth           % 地球图元（Earth/Moon 场景：球面；Synodic 场景：背景散点）
        moon            % 月球图元（Synodic 场景：并入 earth 散点，不再单独使用）
        sats = struct('name', {}, 'line', {}, 'pStart', {}, 'pEnd', {}, ...
            'label', {}, 'color', {}, 'rMax', {})
                        % 已登记卫星：struct 数组（名字、图元句柄、颜色、轨道半径）
        satSerial = 0   % 卫星登记流水号（仅用于缺省命名，删除后不复用）
        labels = struct('h', {}, 'anchor', {}, 'relOff', {})
                        % 文字标注登记表：句柄、锚点、相对偏移（偏移随视图缩放同比变化）
        limListen       % 坐标范围变化监听器（触发文字标注自动重排）
    end

    properties (Constant)
        RE  = 6378.1363;          % 地球赤道半径 (km)
        RM  = 1737.4;             % 月球平均半径 (km)
        REM = 384400;             % 地月平均距离 (km)
        mu  = 0.0121505856313668; % 地月质比 mu = m月/(m地+m月)
                                  % （与 module_const.MOON_EARTH_MASS_RATIO 一致）
        % 卫星轨道配色：固定顺序，不再循环（第 1 个槽位的蓝色留给地球）
        satColors = [ ...
            235 104  52; ...   % 1 橙
             27 175 122; ...   % 2 青
            237 161   0; ...   % 3 黄
            232 123 164]/255;  % 4 品红
    end

    methods
        %>>>>> 构造函数：输入场景名，决定初始化哪个场景
        function obj = viewer(scene)
            % 场景名解析（支持英文与中文别名，缺省为 Earth）
            if nargin < 1 || isempty(scene)
                scene = 'Earth';
            end
            switch lower(scene)
                case {'earth', '3d地球'}
                    obj.scene = 'Earth';
                case {'moon', '3d月球'}
                    obj.scene = 'Moon';
                case {'synodic', '地月会合'}
                    obj.scene = 'Synodic';
                otherwise
                    error('viewer: 未知场景 %s（可选 Earth / Moon / Synodic）', scene);
            end

            % 创建图形窗口和坐标轴
            obj.fig = figure('Color', 'w');
            obj.ax = axes(obj.fig);
            hold(obj.ax, 'on');

            % 按场景绘制背景
            switch obj.scene
                case 'Earth'
                    obj.drawEarth();
                case 'Moon'
                    obj.drawMoon();
                case 'Synodic'
                    obj.drawSynodic();
            end

            % 坐标轴外观：三轴显示比例相同，避免图形变形
            axis(obj.ax, 'equal');                   % 单位长度一致（DataAspectRatio 1:1:1）
            obj.ax.PlotBoxAspectRatio = [1 1 1];     % 绘图盒固定为立方体，防止拉伸填充
            grid(obj.ax, 'on');
            box(obj.ax, 'on');
            obj.ax.GridAlpha = 0.25;
            switch obj.scene
                case 'Earth'
                    obj.ax.Title.String = '卫星轨道视图 (J2000)';
                case 'Moon'
                    obj.ax.Title.String = '卫星轨道视图 (月心惯性系)';
                case 'Synodic'
                    obj.ax.Title.String = '卫星轨道视图 (地月会合系, 归一化)';
            end
            switch obj.scene
                case 'Synodic'
                    obj.ax.XLabel.String = 'X';
                    obj.ax.YLabel.String = 'Y';
                    obj.ax.ZLabel.String = 'Z';
                otherwise
                    obj.ax.XLabel.String = 'X (km)';
                    obj.ax.YLabel.String = 'Y (km)';
                    obj.ax.ZLabel.String = 'Z (km)';
            end
            view(obj.ax, 40, 25);
            rotate3d(obj.ax, 'on');
            obj.refreshLimits();
            % 坐标范围变化时自动重排文字标注（覆盖 zoomTo/zoomOut/工具栏缩放等所有途径）
            obj.limListen = addlistener(obj.ax, {'XLim', 'YLim', 'ZLim'}, ...
                'PostSet', @(src, evt) obj.positionLabels());
        end

        %>>>>> 添加卫星轨道（可多次调用，每次添加一条）
        % 输入：fn 星历文件；name 卫星名字（可选，缺省为 sat1, sat2, ...）
        function touchSat(obj, fn, name)
            D = readEph(fn);

            nCol = size(D, 2);
            if nCol == 8                      % MJD日 + MJD秒 + 位置 + 速度
                r = D(:, 3:5);
            elseif nCol == 7                  % MJD + 位置 + 速度
                r = D(:, 2:4);
            else
                error('viewer: 无法识别的星历格式（%d 列）', nCol);
            end

            % 会合系输入单位为 km、km/s，绘图单位为归一化长度（地月距离）
            if strcmp(obj.scene, 'Synodic')
                r = r / obj.REM;
            end

            obj.satSerial = obj.satSerial + 1;
            if nargin < 3 || isempty(name)
                name = sprintf('sat%d', obj.satSerial);   % 未输入名字时内部定义
            end
            c = obj.pickColor();   % 选取第一个未被占用的配色

            % 轨道线
            hLine = plot3(obj.ax, r(:, 1), r(:, 2), r(:, 3), ...
                'Color', c, 'LineWidth', 1.5);
            % 起点（方块）与终点（圆点）
            p1 = r(1, :);
            pN = r(end, :);
            hStart = plot3(obj.ax, p1(1), p1(2), p1(3), 's', ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', 'none', 'MarkerSize', 9);
            hEnd = plot3(obj.ax, pN(1), pN(2), pN(3), 'o', ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', 'none', 'MarkerSize', 9);

            % 卫星名字标注在起点外侧（偏移随视图缩放同比变化，保持视觉距离不变）
            if strcmp(obj.scene, 'Synodic')
                relD = 0.016;    % 归一化场景（0.02 / 全视图半宽 1.265）
            else
                relD = 500 / (1.1 * max(obj.RE, obj.RM));   % km 场景（保持原 500 km 观感）
            end
            hLabel = text(obj.ax, p1(1), p1(2), p1(3), ...
                name, 'FontSize', 10, 'Color', [0.04 0.04 0.04]);
            obj.registerLabel(hLabel, p1, relD * p1 / norm(p1));

            % 对象内部登记卫星
            obj.sats(end+1) = struct('name', name, 'line', hLine, ...
                'pStart', hStart, 'pEnd', hEnd, 'label', hLabel, ...
                'color', c, 'rMax', max(abs(r(:))));

            % 图例重建：每条卫星一行（仅卫星名）
            obj.updateLegend();

            % 坐标轴范围更新
            obj.refreshLimits();
        end

        %>>>>> 从图中删除卫星
        % 输入：name 卫星名字（与 touchSat 时登记的名字一致）
        function killSat(obj, name)
            idx = find(strcmp({obj.sats.name}, name), 1);
            if isempty(idx)
                warning('viewer: 未找到卫星 %s，未删除', name);
                return;
            end
            % 删除图元：轨道线、起点、终点、名字标注
            delete(obj.sats(idx).line);
            delete(obj.sats(idx).pStart);
            delete(obj.sats(idx).pEnd);
            delete(obj.sats(idx).label);
            % 同步移除文字标注登记
            if ~isempty(obj.labels)
                labIdx = find([obj.labels.h] == obj.sats(idx).label, 1);
                if ~isempty(labIdx)
                    obj.labels(labIdx) = [];
                end
            end
            % 从登记表中移除
            obj.sats(idx) = [];
            % 重建图例，坐标轴范围回缩
            obj.updateLegend();
            obj.refreshLimits();
        end

        %>>>>> 放大观察局部：把视图中心移到目标点并缩小坐标范围
        % 输入：target 目标名（'Earth'/'Moon'/'L1'~'L5'，仅 Synodic 场景）
        %             或 1x3 坐标（Synodic 场景为归一化坐标，其余场景为 km）；
        %       halfRange 窗口半宽（可选；缺省按目标名给定，数值目标缺省 0.1）
        % 说明：三轴使用相同半宽，保持等比例显示；目标点位于视图中心
        function zoomTo(obj, target, halfRange)
            % 解析目标点与缺省窗口半宽
            if isnumeric(target)
                c = target(:)';
                if numel(c) ~= 3
                    error('viewer: 坐标目标须为 1x3 向量');
                end
                h0 = 0.1;
            else
                if ~strcmp(obj.scene, 'Synodic')
                    error('viewer: 命名目标 zoomTo 仅支持 Synodic 场景');
                end
                L = obj.lagrangePoints();
                switch lower(target)
                    case 'earth'
                        c = [-obj.mu, 0, 0];    h0 = 0.05;   % 地球
                    case 'moon'
                        c = [1 - obj.mu, 0, 0]; h0 = 0.02;   % 月球
                    case 'l1'
                        c = L(1, :);            h0 = 0.12;
                    case 'l2'
                        c = L(2, :);            h0 = 0.12;
                    case 'l3'
                        c = L(3, :);            h0 = 0.12;
                    case 'l4'
                        c = L(4, :);            h0 = 0.12;
                    case 'l5'
                        c = L(5, :);            h0 = 0.12;
                    otherwise
                        error('viewer: 未知目标 %s（可选 Earth / Moon / L1~L5）', target);
                end
            end
            if nargin >= 3 && ~isempty(halfRange)
                h0 = halfRange;
            end
            if ~(h0 > 0)
                error('viewer: 窗口半宽须为正数');
            end
            % 以目标点为中心重建坐标范围（三轴同宽，保持等比例）
            obj.ax.XLim = c(1) + [-h0, h0];
            obj.ax.YLim = c(2) + [-h0, h0];
            obj.ax.ZLim = c(3) + [-h0, h0];
            obj.positionLabels();
        end

        %>>>>> 恢复全景：坐标范围覆盖背景与所有轨道
        function zoomOut(obj)
            obj.refreshLimits();
        end

        %>>>>> 析构：清理坐标范围监听器
        function delete(obj)
            delete(obj.limListen);
        end
    end

    methods (Access = private)
        %>>>>> 场景一：3D 地球（worldMap1.jpg 贴图）
        %      贴图来源：NASA Blue Marble（NASA Earth Observatory / GSFC），详见 THIRD_PARTY.md
        function drawEarth(obj)
            mapFile = fullfile(fileparts(mfilename('fullpath')), 'worldMap1.jpg');
            [x, y, z] = sphere(50);
            if exist(mapFile, 'file') == 2
                obj.earth = surf(obj.ax, x*obj.RE, y*obj.RE, z*obj.RE, ...
                    'EdgeColor', 'none', 'FaceColor', 'texturemap', ...
                    'CData', imread(mapFile));
                %>>>>> 旋转和翻转，调整贴图，使它与经纬度对应
                %      150° 为经验值，若 worldMap1.jpg 经度对不齐可在此微调
                set(obj.earth, 'ZData', -get(obj.earth, 'ZData')); % Z方向翻转，使北半球朝上
                rotate(obj.earth, [0 0 1], 150);                    % 沿Z轴旋转，对齐经度
            else
                % 找不到贴图时退化为纯色球体
                obj.earth = surf(obj.ax, x*obj.RE, y*obj.RE, z*obj.RE, ...
                    'EdgeColor', 'none', 'FaceColor', [158 197 244]/255);
                camlight(obj.ax, 'headlight');
                lighting(obj.ax, 'gouraud');
                warning('viewer: 未找到贴图 %s，改用纯色地球', mapFile);
            end
            obj.drawEquator(obj.RE);
        end

        %>>>>> 场景二：3D 月球（灰色球体 + 光照）
        function drawMoon(obj)
            [x, y, z] = sphere(50);
            obj.earth = surf(obj.ax, x*obj.RM, y*obj.RM, z*obj.RM, ...
                'EdgeColor', 'none', 'FaceColor', [0.55 0.55 0.55]);
            camlight(obj.ax, 'headlight');
            lighting(obj.ax, 'gouraud');
            obj.drawEquator(obj.RM);
        end

        %>>>>> 场景三：地月会合坐标系（归一化单位，旋转系）
        % 长度单位 = 地月距离；原点 = 地月质心；
        % 地球 [-mu, 0, 0]，月球 [1-mu, 0, 0]，静止不动；
        % 地球、月球、L1~L5 均为等大小散点（淡蓝 / 灰 / 深灰）
        function drawSynodic(obj)
            L = obj.lagrangePoints();
            % 地球、月球与 L1~L5：一个散点图，等大小
            P = [-obj.mu, 0, 0; ...    % 地球
                 1 - obj.mu, 0, 0; ... % 月球
                 L];                   % L1~L5
            C = [0.529 0.808 0.922; ...              % 地球：淡蓝
                 0.550 0.550 0.550; ...              % 月球：灰
                 repmat([0.250 0.250 0.250], 5, 1)]; % L1~L5：深灰
            S = [45, 45, 18, 18, 18, 18, 18];        % 地球/月球：90 的 1/2；L1~L5：90 的 1/5
            obj.earth = scatter3(obj.ax, P(:, 1), P(:, 2), P(:, 3), S, C, ...
                'filled', 'MarkerEdgeColor', 'none');
            % 文字说明（锚点为天体位置，偏移随视图缩放同比变化）
            hT = text(obj.ax, -obj.mu, 0, 0, 'Earth', ...
                'FontSize', 11, 'Color', [0.04 0.04 0.04], ...
                'HorizontalAlignment', 'center');
            obj.registerLabel(hT, [-obj.mu, 0, 0], [0.063, 0.063, 0.063]);
            hT = text(obj.ax, 1 - obj.mu, 0, 0, 'Moon', ...
                'FontSize', 11, 'Color', [0.04 0.04 0.04], ...
                'HorizontalAlignment', 'center');
            obj.registerLabel(hT, [1 - obj.mu, 0, 0], [0.063, 0.063, 0.063]);
            % 地月连线（灰色虚线，标示会合系 +X 方向）
            plot3(obj.ax, [-obj.mu, 1 - obj.mu], [0, 0], [0, 0], '--', ...
                'Color', [137 135 129]/255, 'LineWidth', 0.5);
            % 拉格朗日点 L1~L5 文字说明（偏移随视图缩放同比变化）
            labL = {'L1', 'L2', 'L3', 'L4', 'L5'};
            relL = [0.047, 0.071, 0.055];   % 原固定偏移 (0.06, 0.09, 0.07) / 全视图半宽 1.265
            for k = 1:5
                hT = text(obj.ax, L(k,1), L(k,2), L(k,3), labL{k}, ...
                    'FontSize', 8, 'Color', [0.04 0.04 0.04], ...
                    'HorizontalAlignment', 'center');   % L1~L5 文字：原大小 1/2
                obj.registerLabel(hT, L(k, :), relL);
            end
            % L4-L5 连线（灰色虚线，便于观察两点相对位置）
            plot3(obj.ax, L(4:5, 1), L(4:5, 2), L(4:5, 3), '--', ...
                'Color', [137 135 129]/255, 'LineWidth', 0.5);
        end

        %>>>>> 拉格朗日点 L1~L5 位置（归一化会合系；公式与 module_CRTBP 一致）
        function L = lagrangePoints(obj)
            mu0 = obj.mu;
            f = 1/3;
            g = (mu0*f)^f;
            xi1 = g * (1 - f*mu0*g - f/3*g^2);
            xi2 = g * (1 + f*mu0*g - f/3*g^2);
            nu = 7/12*mu0;
            xi3 = 1 - nu*(1 + 23/84*nu^2 + 23/84*nu^3 + ...
                          761/2352*nu^4 + 3163/7056*nu^5 + 30703/49392*nu^6);
            L = [1 - mu0 - xi1, 0, 0; ...        % L1
                 1 - mu0 + xi2, 0, 0; ...        % L2
                 -mu0 - xi3, 0, 0; ...           % L3
                 0.5 - mu0,  sqrt(3)/2, 0; ...   % L4
                 0.5 - mu0, -sqrt(3)/2, 0];      % L5
        end

        %>>>>> 赤道参考圈（灰色虚线）
        function drawEquator(obj, R)
            th = linspace(0, 2*pi, 361);
            plot3(obj.ax, R*cos(th), R*sin(th), zeros(size(th)), ...
                '--', 'Color', [137 135 129]/255, 'LineWidth', 0.5);
        end

        %>>>>> 重建图例：每条卫星一行（仅卫星名）
        function updateLegend(obj)
            hLeg = [];
            labLeg = {};
            for k = 1:numel(obj.sats)
                hLeg = [hLeg, obj.sats(k).line]; %#ok<AGROW>
                labLeg = [labLeg, obj.sats(k).name]; %#ok<AGROW>
            end
            legend(obj.ax, hLeg, labLeg, 'Location', 'northeast');
        end

        %>>>>> 选取第一个未被占用的配色；全部占用时用灰色
        function c = pickColor(obj)
            used = false(size(obj.satColors, 1), 1);
            for k = 1:numel(obj.sats)
                for m = 1:size(obj.satColors, 1)
                    if isequal(obj.sats(k).color, obj.satColors(m, :))
                        used(m) = true;
                        break;
                    end
                end
            end
            m = find(~used, 1);
            if isempty(m)
                c = [137 135 129]/255;   % 颜色耗尽：灰色（不再生成新颜色）
            else
                c = obj.satColors(m, :);
            end
        end

        %>>>>> 更新坐标轴范围：覆盖背景天体与所有已登记轨道
        function refreshLimits(obj)
            switch obj.scene
                case 'Earth'
                    maxR = obj.RE;
                case 'Moon'
                    maxR = obj.RM;
                case 'Synodic'
                    maxR = 1.15;   % 覆盖地球、月球与 L1~L5（L2 附近轨道加入后自动扩大）
            end
            for k = 1:numel(obj.sats)
                maxR = max(maxR, obj.sats(k).rMax);
            end
            lim = 1.1 * maxR;
            obj.ax.XLim = [-lim, lim];
            obj.ax.YLim = [-lim, lim];
            obj.ax.ZLim = [-lim, lim];
            obj.positionLabels();
        end

        %>>>>> 登记文字标注：句柄 + 锚点 + 相对偏移
        % relOff 为相对当前视图半宽的比例，保证放大/缩小后视觉距离不变
        function registerLabel(obj, hText, anchor, relOff)
            obj.labels(end+1) = struct('h', hText, 'anchor', anchor, ...
                'relOff', relOff);
        end

        %>>>>> 按当前坐标范围重排文字标注（偏移 = relOff × 视图半宽）
        function positionLabels(obj)
            h = max([diff(obj.ax.XLim), diff(obj.ax.YLim), diff(obj.ax.ZLim)]) / 2;
            for k = 1:numel(obj.labels)
                set(obj.labels(k).h, 'Position', ...
                    obj.labels(k).anchor + obj.labels(k).relOff * h);
            end
        end
    end
end

%>>>>> 读取星历（跳过 # 注释行，按行解析为数值）
function D = readEph(fn)
    if exist(fn, 'file') ~= 2
        error('viewer: 找不到星历文件 %s', fn);
    end
    fid = fopen(fn, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = lines{1};
    D = cellfun(@(s) sscanf(s, '%f')', lines, 'UniformOutput', false);
    D = cell2mat(D(~cellfun(@isempty, D)));
end
