function [bytes, source] = available_memory_bytes()
%AVAILABLE_MEMORY_BYTES  Physical memory currently available to new allocations, in bytes.
%
%   [bytes, source] = AVAILABLE_MEMORY_BYTES()
%
%   ADDED 2026-10-27 for the memory check of run_reliability_draws. Returns NaN (source
%   'unknown') when the platform gives no answer; callers then skip their check.
%     Windows  memory(): PhysicalMemory.Available (MATLAB, and Octave where memory() exists)
%     Linux    /proc/meminfo: MemAvailable
%     macOS    vm_stat: (free + inactive + speculative pages) x page size
%   The figure is physical memory only: swap is not counted, so a run that passes the check
%   does not have to page.

bytes = NaN;  source = 'unknown';
try
    if ispc && exist('memory') > 0                                          %#ok<EXIST>
        [~, s] = memory;
        bytes = double(s.PhysicalMemory.Available);  source = 'memory()';
    elseif exist('/proc/meminfo', 'file')
        txt = fileread('/proc/meminfo');
        tok = regexp(txt, 'MemAvailable:\s*(\d+)\s*kB', 'tokens', 'once');
        if ~isempty(tok), bytes = str2double(tok{1}) * 1024;  source = '/proc/meminfo'; end
    elseif ismac
        [status, out] = system('vm_stat');
        if status == 0
            page = regexp(out, 'page size of (\d+) bytes', 'tokens', 'once');
            pages = 0;
            for name = {'Pages free', 'Pages inactive', 'Pages speculative'}
                tok = regexp(out, [name{1} ':\s*(\d+)'], 'tokens', 'once');
                if ~isempty(tok), pages = pages + str2double(tok{1}); end
            end
            if ~isempty(page), bytes = pages * str2double(page{1});  source = 'vm_stat'; end
        end
    end
catch
    bytes = NaN;  source = 'unknown';
end
end
