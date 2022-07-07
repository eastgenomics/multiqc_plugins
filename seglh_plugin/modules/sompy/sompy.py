
#!/usr/bin/env python

""" MultiQC som.py plugin module """

from __future__ import print_function
from collections import OrderedDict, defaultdict
import logging
import os
import re
from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')


def autocast(x):
    '''automatically typecasts numerical values to int or float'''
    if re.match('^(\d+)\.(\d+)$', x):
        return float(x)
    elif re.match('^(\d+)$', x):
        return int(x)
    else:
        return x


class MultiqcModule(BaseMultiqcModule):
    # configures manual metric->group mappings
    sompy_groups = {
        'records':'Combined Benchmark',
        'SNVs': 'SNV Benchmark',
        'indels':'Indel Benchmark'
    } 

    def __init__(self):
        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name = 'som.py metrics',
            target = "sompy",
            anchor = 'sompy',
            href = 'https://github.com/moka-guys/multiqc_plugins',
            info = "som.py variant call benchmarking."
        )

        # Find and load any input files for this module
        self.sompy_data = defaultdict(dict)
        self.source_files = dict()
        for f in self.find_log_files('sompy'):
            self.parse_file(f)
            self.add_data_source(
                s_name=f['s_name'],
                source=os.path.join(f['root'],f['fn']),
                module="sompy",
                section="sompy-bysample",
            )

        # Filter out samples matching ignored sample names
        self.sompy_data = self.ignore_samples(self.sompy_data)

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.sompy_data) == 0:
            log.debug("Could not find any som.py reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.sompy_data)))

        # Write parsed report data to a file
        combined_data = dict()
        for sample in self.sompy_data:
            combined_data[sample] = dict()
            for group in self.sompy_data[sample]:
                for metric in self.sompy_data[sample][group]:
                    combined_data[sample][f'{group}_{metric}'] = self.sompy_data[sample][group][metric]
        self.write_data_file(combined_data, 'multiqc_sompy')

        # create the result table
        for group in sorted(self.sompy_groups.keys()):
            self.add_section(
                name=self.sompy_groups[group],
                anchor="sompy-bysample",
                description="",
                plot=self.sample_stats_table(group),
            )

    def sample_stats_table(self, group):
        '''
        create a table with the sample statistics
        '''
        h = OrderedDict()
        h["unk"] = {
            "title": "Not assessed calls",
            "description": "Number of non-assessed query calls",
            "min": 0,
            "max": 1,
            "hidden": True,
            "format": "{:.4f}",
        }

        h["total.truth"] = {
            "title": "Truth: Total",
            "description": "Total number of truth variants",
            "format": None,
            "hidden": True,
        }

        h["total.query"] = {
            "title": "Query: Total",
            "description": "Total number of query calls",
            "format": None,
            "hidden": True,
        }

        h["tp"] = {
            "title": "True Positives",
            "description": "Number of true-positive calls",
            "suffix": " variants",
            "scale": "Reds",
            "format": None,
        }

        h["fn"] = {
            "title": "False Negatives",
            "description": "Calls in truth without matching query call",
            "suffix": " variants",
            "scale": "Reds",
            "format": None,
        }
        h["fp"] = {
            "title": "False Positives",
            "description": "Number of false-positive calls",
            "format": None,
            "scale": "Reds",
            "hidden": True,
        }
        h["unk"] = {
            "title": "Unknown",
            "description": "Number of calls outside the confident regions",
            "format": None,
            "hidden": True,
        }
        h["recall"] = {  # string must match headers in the input file
            "title": "Recall",  # whatever string to be displayed in the html report table
            "description": "Recall for truth variant representation = TRUTH.TP / (TRUTH.TP + TRUTH.FN)",
            "min": 0,
            "max": 1,
            "cond_formatting_rules": {
                "verygreen": [
                    {"gte": 0.99},
                ],
                "green": [
                    {"lt": 0.99},
                    {"gt": 0.98}
                ],
                "amber": [
                    {"lt": 0.98},
                    {"gt": 0.90}
                ],
                "red": [
                    {"lt": 0.90}
                ],
            },
            "cond_formatting_colours": [
                {"red": "#D2222D"},
                {"amber": "#FFBF00"},
                {"green": "#238823"},
                {"verygreen": "#007000"},
            ],
            "format": "{:.4f}",
        }
        h["precision"] = {
            "title": "Precision",
            "description": "Precision of query variants = QUERY.TP / (QUERY.TP + QUERY.FP)",
            "min": 0,
            "max": 1,
            "format": "{:.4f}",
            "hidden": True,
        }

        # Table config
        table_config = {
            "namespace": "sompy",
            "id": f"sompy-sample-stats-table-{group}",
            "table_title": f"sompy Sample Statistics ({group})",
            "no_beeswarm": True,
        }
        group_data = dict()
        for sample in self.sompy_data:
            group_data[sample] = self.sompy_data[sample][group]

        return table.plot(group_data, h, table_config)


    def parse_file(self, f):
        '''Parses the Metrics output file
        CSV file with header, one row per variant type

        input:
            f: file handle
        output:
            None
        '''
        header = []
        group, sample_names = '', []
        for line in f['f'].splitlines():
            if line.endswith('sompyversion,sompycmd'):
                # header
                header = line.rstrip().split(',') 
            elif line.startswith('#') or len(line) == 0 or re.match(r'^\s+$',line):
                # comment or empty line (reset)
                continue
            else:
                # data line, extract and typecast numbers
                fields = list(map(autocast, line.rstrip().split(',')))
                # add sample data dictionaries
                group = fields[1]
                data = dict(zip(header[2:], fields[2:]))
                # get sample identifier from output name
                sample_name = None
                # from the sompy file, select the column sompy_cmd and
                # matches for -o {output file}
                outfile = re.search(r'.*-o\s*(\S+)', data['sompycmd'])
                if outfile:
                    sample_name = os.path.basename(outfile.group(1))
                    # add to data dictionary
                    self.sompy_data[sample_name][group] = data
                else:
                    log.debug("Could not find string -o in sompy stats.csv file")
                    raise UserWarning
