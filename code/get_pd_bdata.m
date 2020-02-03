% Gets protocol_data using bdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% CHUCK RATS, FROM DATA_SETTINGS_TRAWL.JL

rats = {'select protocol_data from sessions where ratname="K238" and (sessiondate>="2017-09-11" and sessiondate<="2017-10-02") order by sessiondate',
        'select protocol_data from sessions where ratname="K241" and (sessiondate>="2017-11-27" and sessiondate<="2017-12-13") order by sessiondate',
        'select protocol_data from sessions where ratname="K265" and (sessiondate>="2019-06-12" and sessiondate<="2019-07-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K283" and (sessiondate>="2019-07-23" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K284" and (sessiondate>="2019-07-23" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K285" and (sessiondate>="2019-07-23" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K289" and (sessiondate>="2019-07-23" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K294" and (sessiondate>="2018-10-13" and sessiondate<="2018-10-18") order by sessiondate',
        'select protocol_data from sessions where ratname="K296" and (sessiondate>="2018-11-10" and sessiondate<="2020-01-16") order by sessiondate',
        'select protocol_data from sessions where ratname="K299" and (sessiondate>="2018-10-16" and sessiondate<="2018-10-18") order by sessiondate',
        'select protocol_data from sessions where ratname="K310" and (sessiondate>="2019-07-23" and sessiondate<="2019-08-13") order by sessiondate',
        'select protocol_data from sessions where ratname="K311" and (sessiondate>="2019-07-23" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K313" and (sessiondate>="2019-07-23" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K314" and (sessiondate>="2019-07-23" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K316" and (sessiondate>="2019-07-23" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K317" and (sessiondate>="2019-08-18" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K319" and (sessiondate>="2019-07-23" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K322" and (sessiondate>="2019-07-23" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K323" and (sessiondate>="2019-12-19" and sessiondate<="2020-01-24") order by sessiondate',
        'select protocol_data from sessions where ratname="K324" and (sessiondate>="2019-12-02" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K328" and (sessiondate>="2019-09-05" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K329" and (sessiondate>="2019-11-08" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K330" and (sessiondate>="2019-11-03" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K331" and (sessiondate>="2019-10-30" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K332" and (sessiondate>="2019-12-04" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K335" and (sessiondate>="2019-10-18" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K336" and (sessiondate>="2019-10-03" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K338" and (sessiondate>="2019-12-23" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K339" and (sessiondate>="2019-09-25" and sessiondate<="2020-02-01") order by sessiondate',
        'select protocol_data from sessions where ratname="K340" and (sessiondate>="2020-01-17" and sessiondate<="2020-01-31") order by sessiondate',
        'select protocol_data from sessions where ratname="K341" and (sessiondate>="2019-11-11" and sessiondate<="2020-01-31") order by sessiondate'};

tasktype = {'classic',
            'classic',
            'classic',
            'classic',
            'classic',
            'classic',
            'classic',
            'frequency',
            ['frequency' 'classic'],
            'frequency',
            ['frequency' 'classic'],
            'frequency',
            'frequency',
            'frequency',
            'frequency',
            'frequency',
            'classic',
            'classic',
            'classic',
            'classic',
            'classic',
            'classic',
            'classic',
            'classic',
            'frequency',
            'classic',
            'frequency',
            'classic',
            'classic',
            'frequency',
            'frequency'};

rates = {
    [40.0, 1.0]
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0, 20.0],
    [40.0, 20.0],
    [40.0],
    [40.0, 20.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0, 20.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [40.0],
    [20.0],
    [40.0],
    [20.0],
    [40.0],
    [40.0],
    [20.0],
    [40.0],
    [20.0]
};

base_freq = {
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [[6500.0, 14200.0]],
    [[6500.0, 14200.0], 2000.0],
    [[6500.0, 14200.0]],
    [[6500.0, 14200.0], 2000.0, [4500.0, 10400.0]],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [2000.0],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [2000.0],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [2000.0],
    [2000.0],
    [[6500.0, 14200.0], [4500.0, 10400.0]],
    [[6500.0, 14200.0], [4500.0, 10400.0]]
};

cfg.EXPORTPATH_DATA = './';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for irat = 1 : length(rats)
    pds = bdata(rats{irat});
    if length(pds{1}.hits) == length(pds{1}.bupsdata)
        ratname = split(rats{irat}, '"');
        ratname = ratname{2};
        fprintf('Running! %s\n', ratname)
        % ratdata = parse_pd_sessions(pds);
        hh = []; bt = [];
        vh = []; bd = [];
        g  = []; sd = [];
        gr = []; b = [];

        for isess = 1 : length(pds)
            hh = [hh ; pds{isess}.hits];
            bt = [bt ; pds{isess}.bupsdata];
            vh = [vh ; pds{isess}.violations];
            bd = [bd ; pds{isess}.n_right - pds{isess}.n_left];
            sd = [sd ; pds{isess}.samples];

            g_pd = [];
            for itrl = 1 : length(pds{isess}.hits)
                g_pd = [g_pd ; pds{isess}.bupsdata{itrl}.gamma];
            end
            g = [g ; g_pd];

            gr_pd = [];
            for itrl = 1 : length(pds{isess}.hits)
                if ~isnan(pds{isess}.hits(itrl))
                    gr_pd = [gr_pd ; (pds{isess}.hits(itrl) && pds{isess}.sides(itrl) == 'r') ||...
                                     (pds{isess}.hits(itrl) == 0 && pds{isess}.sides(itrl) == 'l')];
                else
                    gr_pd = [gr_pd ; NaN];
                end
            end
            gr = [gr ; gr_pd];

            % b_pd = [];
            for itrl = 1 : length(pds{isess}.hits)
                b(end+1).left = pds{isess}.bupsdata{itrl}.left;
                b(end).right  = pds{isess}.bupsdata{itrl}.right;
            end
            % b(end+1).left = b_pd.left;
            % b(end).right  = b_pd.right;
        end

        ratdata = [];
        ratdata.hh = hh;
        ratdata.bt = bt;
        ratdata.vh = vh;
        ratdata.bd = bd;
        ratdata.g  = g;
        ratdata.sd = sd;
        ratdata.gr = gr;
        ratdata.b  = b;

        % Now let's generate the parsed field, removing all the NaNs
        goodtrls = ~isnan(ratdata.hh);
        ratdata.parsed = struct();
        ratdata.parsed.hh = ratdata.hh(goodtrls);
        ratdata.parsed.bd = ratdata.bd(goodtrls);
        ratdata.parsed.vh = ratdata.vh(goodtrls);
        ratdata.parsed.bd = ratdata.bd(goodtrls);
        ratdata.parsed.g  = ratdata.g(goodtrls);
        ratdata.parsed.sd = ratdata.sd(goodtrls);
        ratdata.parsed.gr = ratdata.gr(goodtrls);
        ratdata.parsed.pd = ratdata.sd(goodtrls);
        ratdata.parsed.b  = ratdata.b(goodtrls);

        save([cfg.EXPORTPATH_DATA ratname '_' tasktype{irat} '_' num2str(rates{irat}) 'Hz_' num2str(length(base_freq{irat})) '.mat'], 'ratdata');
    end
end

% """
% ratdata = 
%   struct with fields:
%         hh: [1×50036 double]
%         bt: {1×50036 cell}
%         vh: [1×50036 double]
%         vc: []
%         pd: [1×83 double]
%         sd: [1×50036 double]    | pds{i}.bupsdata{j}.real_T
%         sh: [1×50036 double]
%         gr: [1×50036 double]
%       days: {}
%        sid: []
%     parsed: [1×1 struct]
%      psych: {[1×1 struct]}
% """
function ratdata = parse_pd_sessions(pds)
    hh = []; bt = [];
    vh = []; bd = [];
    g  = []; sd = [];
    gr = [];

    for isess = 1 : length(pds)
        hh = [hh ; pds{isess}.hits];
        bt = [bt ; pds{isess}.bupsdata];
        vh = [vh ; pds{isess}.violations];
        bd = [bd ; pds{isess}.n_right - pds{isess}.n_left];
        sd = [sd ; pds{isess}.samples];

        g_pd = [];
        for itrl = 1 : length(pds{isess}.hits)
            g_pd = [g_pd ; pds{isess}.bupsdata{itrl}.gamma];
        end
        g = [g ; g_pd];

        gr_pd = [];
        for itrl = 1 : length(pds{1}.hits)
            if ~isnan(pds{1}.hits(itrl))
                gr_pd = [gr_pd ; (pds{1}.hits(itrl) && pds{1}.sides(itrl) == 'r') || (pds{1}.hits(itrl) == 0 && pds{1}.sides(itrl) == 'l')];
            else
                gr_pd = [gr_pd ; NaN];
            end
        end
        gr = [gr ; gr_pd];

        % sds = [];
        % for itrl = 1 : length(pds{isess}.hits)
        %     sds = [sds ; pds{isess}.bupsdata{itrl}.real_T];
        % end
        % sd = [sd ; sds];
    end

    % disp(size(hh));
    % disp(size(bt));
    % disp(size(vh));
    % disp(size(sd));
    % disp(size(gr));

    ratdata = [];
    ratdata.hh = hh;
    ratdata.bt = bt;
    ratdata.vh = vh;
    ratdata.bd = bd;
    ratdata.g  = g;
    ratdata.sd = sd;
    ratdata.gr = gr;

    % ratdata = struct('hh', hh, 'bt', bt, 'vh', vh, 'vc', [], 'pd', sd,...
    %                  'sd', sd, 'sh', [], 'gr', gr, 'days', {}, 'sid', [],...
    %                  'parsed', [], 'psych', {});

    % ratdata.hh = hh;
    % ratdata.bt = bt;
    % ratdata.vh = vh;
    % ratdata.bd = bd;
    % ratdata.g  = g ;
    % ratdata.sd = sd;
    % ratdata.gr = gr;

    % Now let's generate the parsed field, removing all the NaNs
    goodtrls = ~isnan(ratdata.hh);
    ratdata.parsed = struct();
    ratdata.parsed.hh = ratdata.hh(goodtrls);
    ratdata.parsed.bd = ratdata.bd(goodtrls);
    ratdata.parsed.vh = ratdata.vh(goodtrls);
    ratdata.parsed.bd = ratdata.bd(goodtrls);
    ratdata.parsed.g  = ratdata.g(goodtrls);
    ratdata.parsed.sd = ratdata.sd(goodtrls);
    ratdata.parsed.gr = ratdata.gr(goodtrls);
    ratdata.parsed.pd = ratdata.sd(goodtrls);
end


