import orchestration.io as cpio
import lib.osutils as osutil
from lib.fitsutil import createFitsFromDB
from tools.gavo.dao import CatalogDB
import pyfits
import ephem
import os
import glob
import re
import orchestration
import tools.base_component
import lib.config
import lib.log
import pyfits
import random

mwfitting_package = "mwfitting 1.21+0"

mw_fitting_map = [
        mwfitting_package,
        "matplotlib 1.5.2+0",
        "pyfits 3.4+0",
        "trilegal 1.6+0",
        "healpy 1.9.1+1",
        "astropy 1.2.1+0",
        "Python 3.4.3+0"
]

mw_fitting_multi_make_hess = [
        mwfitting_package,
        "trilegal 1.6+0",
        "numpy 1.7.1+1",
        "pyfits 3.3+1",
        "pprocess 0.5.1+0",             # bad setuptools -> bad python
        "setuptools 14.0+0",
        "python 2.7.3+2"
]

mw_fitting_modelfit = [
        mwfitting_package,
        "trilegal 1.6+0",
        "emcee 2.1.0+0",
        "scipy 0.11.0+2",
        "matplotlib 2.0.2+0",
        "numpy 1.7.1+1",
        "corner 2.0.1+0",
#        "astropy 1.0.4+1",             # lot of uncompatible deps, better to use old pyfits
        "pyfits 3.3+1",
        "pprocess 0.5.1+0",             # bad setuptools -> bad python
        "setuptools 14.0+0",
        "python 2.7.3+2"
]

mw_fitting_multi_make_mock = [
        mwfitting_package,
        "trilegal 1.6+0",
        "numpy 1.7.1+1",
        "pyfits 3.3+1",
        "python 2.7.3+2"
]

initial_python_path = os.path.join(os.getenv("ASTROSOFT_ROOT"), "pypeline")

class MwFitting(tools.base_component.BaseComponent):
    def run(self, debug=True):
        conf = cpio.ComponentConfig()

        pio = cpio.ProvenanceIO()
        vac = pio.product_by_class.get("vac_ga", None)

        if not vac:
            vac = pio.product_by_class.get("vac_generic", None)
            if not vac:
                raise IOError("No selected products with vac_ga or vac_generic classes.")
 
        # footprint_map = self.create_footprint_map_fits(vac)
        # vac_ga = self.create_vac_ga_fits()

        footprint_map = str()
        vac_ga = str()

        config_file = self.create_config(self.conf, footprint_map, vac_ga)

        self.execute_map()
        self.execute_makehess()
        self.execute_emcee()
        self.execute_multi_make_mock()
        self.execute_dist_stars()
        self.create_product_log()

    def create_config(self, config_object, footprint_map, vac_ga):
        """ create config file """
        self.info("\n\n\n= Creating config file =")

        pathreal = os.getcwd()
        self.debug("##### path real ##### %s" % pathreal)

        vac_ga = "%s/%s" % (pathreal, vac_ga)
        footprint_map = "%s/%s" % (pathreal, footprint_map)
        schlegeldir = "/archive/external_catalogs/SDSS3/schlegel/schlegel"
        datadir = "%s/data" % pathreal
        plotdir = "%s/plots" % pathreal

        if not os.path.exists(datadir):
            os.makedirs(datadir)

        if not os.path.exists(plotdir):
            os.makedirs(plotdir)

        _apply_fitparams = config_object.getCheckedDictByType(
            'APPLY_FITPARAMS').keys()
        apply_fitparams = [x[1:] for x in _apply_fitparams]

        subsections = config_object.doc.xpath('.//section/subsection')

        lines = []

        lines.append("[EXTINCTION]\n")
        lines.append("schlegeldir = %s\n" % schlegeldir)

        for subsection in subsections:
            subsection_id = subsection.attributes[(None, 'id')].nodeValue
            lines.append("[%s]\n" % subsection_id)

            if (subsection_id == "INPUTDATA"):
                #lines.append("vac_ga = %s /archive/external_catalogs/SDSS3/DR14/SDSS_DR14_stars_allgri.fits\n" % vac_ga)
                #lines.append("footprint_map = %s /archive/external_catalogs/SDSS3/DR14/SDSS_DR14_ftp_4096.fits\n" % footprint_map)
                lines.append("vac_ga = /home/adriano.pieres/DES_DR1/DES_DR1_stars.fits /archive/external_catalogs/SDSS3/DR14/SDSS_DR14_stars_allgri.fits\n")
                lines.append("footprint_map = /home/adriano.pieres/DES_DR1/DES_DR1_ftp_512.fits /home/adriano.pieres/SDSS_DR14/SDSS_DR14_ftp_512.fits\n")

            if (subsection_id == "CONFIGHESS"):
                lines.append("datadir = %s\n" % datadir)

            if (subsection_id == "OPTIMISE"):
                lines.append("plotdir = %s\n" % plotdir)

            for scalar in subsection.xpath('./scalar'):
                key = scalar.attributes[(None, 'id')].nodeValue
                value = scalar.attributes[(None, 'value')].nodeValue
                if (subsection_id == "FITPARAMS"):
                    if key in apply_fitparams:
                        lines.append(
                            "%s = %s\n" % (key, value.replace(';', ' ')))
                elif (subsection_id == "SURVEY") and (key == 'survey_name'):
                        lines.append("%s = %s\n" % (key, value.replace(';', ' ')))
                elif (subsection_id == "CONFIGHESS") and (
                                key == 'mag_range' or key == 'color_range' or key =='maxfields' or key == 'mag_cols'):
                        lines.append(
                                "%s = %s\n" % (key, value.replace(';', ' ')))
                else:
                    lines.append("%s = %s\n" % (key, value))

            lines.append("\n")

        config_ini = "mwf.ini"
        _file = open(config_ini, 'w')
        _file.writelines(lines)
        _file.close()

        return config_ini

    def create_vac_ga_fits(self):
        """ create vac_ga fits """
        self.info("\n\n\n= Creating vac ga fits =")

        pio = cpio.ProvenanceIO()

        vac = pio.product_by_class.get("vac_ga", None)

        if not vac:
            vac = pio.product_by_class.get("vac_generic", None)
            if not vac:
                raise IOError("No selected products with vac_ga or vac_generic classes.")

        schema = vac.get('schema')
        table = vac.get('table')

        vac_product = "%s.%s" % (schema, table)

        self.info("vac ga table: %s" % vac_product)

        vac_file = self.create_fits_file(vac_product)

        if not os.path.isfile(vac_file):
            raise IOError('%s file does not exist' % vac_file)

        hdulist = pyfits.open(vac_file)

        tb = hdulist[1].data
        columns = hdulist[1].columns

        cols = []

        col_l = []
        col_b = []

        for col in columns:
            if col.name == 'l':
                col_l = None

            if col.name == 'b':
                col_b = None

            cols.append(
                pyfits.Column(
                    name=col.name,
                    format=col.format,
                    array=tb.field(col.name)
                )
            )

        if col_l is not None or col_b is not None:
            for item in tb:
                if col_l is not None:
                    col_l.append(ephem.Galactic(
                        ephem.Equatorial(item['ra'], item['dec'])).long)

                if col_b is not None:
                    col_b.append(ephem.Galactic(
                        ephem.Equatorial(item['ra'], item['dec'])).lat)

            if col_b is not None:
                col_b = pyfits.Column(name='b', format='D', array=col_b)

            if col_l is not None:
                col_l = pyfits.Column(name='l', format='D', array=col_l)

            if col_l is not None:
                cols = pyfits.ColDefs(cols).add_col(col_l)

            if col_b is not None:
                cols = pyfits.ColDefs(cols).add_col(col_b)

            hdu = pyfits.new_table(cols)
            hdu.writeto(vac_file, clobber=True)

        self.info("--- Partial run time: %s" % self.get_partial_time())
        self.reset_partial_time()
        self.info("--- Comulative run time: %s\n\n" % self.get_total_time())

        return vac_file

    def create_footprint_map_fits(self, vac):
        """ create footprint_map fits """
        self.info("\n\n\n= Creating footprint map fits =")

        vacpid = vac['process_id']

        # get vac path
        vacpath = os.path.join(vac['process_path'], 'pipe_output.xml')

        if not os.path.exists(vacpath):
            raise IOError((
                "vac pipe_output.xml not "
                "exists (ID: %s - PATH: %s)"
            ) % (vacpid, vacpath))

        pout = cpio.ProcessOutputIO(vacpath)

        footprint_table = pout.get_maps_by_class("footprint_map")

        if len(footprint_table) < 1:
            message = (
                "footprint map was not found. " +
                "(pipe output path: %s)" % vacpath)
            self.error(message)
            raise Exception(message)

        if len(footprint_table) > 1:
            self.warn(
                "were found more than a footprint map, " +
                "and is being used the last declared map. " +
                "Footprint maps: %s)" % ", ".join(footprint_table))

        footprint_table = footprint_table.pop()

        self.info("footprint table: %s" % footprint_table)

        footprint_map_file = self.create_fits_file(footprint_table)

        if not os.path.isfile(footprint_map_file):
            raise IOError('%s file does not exist' % footprint_map_file)

        self.info("--- Partial run time: %s" % self.get_partial_time())
        self.reset_partial_time()
        self.info("--- Comulative run time: %s\n\n" % self.get_total_time())

        return footprint_map_file

    def create_fits_file(self, tablename):
        """ create fits file by tablename """

        db = CatalogDB()

        schema, table = tablename.split(".")

        if not db.table_exists(schema, table):
            message = "%s was not found in gavo database. " % tablename
            self.error(message)
            raise Exception(message)

        query = (
            "SELECT * FROM %s " % tablename
        )
        query_result = db.executeIndependentQuery(query)

        fits_file = tablename + ".fits"

        self.info("Creating %s based in: %s" % (fits_file, query))

        createFitsFromDB(
            tablename,
            fits_file,
            query_result,
            db,
            mode="overwrite",
            product_info=None
        )

        db.commit()

        self.info("%s created" % fits_file)

        return fits_file

    def execute_map(self):
        """ execute map """
        self.info("\n\n\n= Running map =")

        cmd = (
            "ulimit -s unlimited && "
            #"python3 -u ${MWFITTING_DIR}/map.py &> map.log"
            "python3 -u /home/adriano.pieres/mwf_git/mwfitting/map.py &> map.log"
        )

        try:
            osutil.OSUtils().run_eups_command(
                cmd,
                mw_fitting_map,
                stderrAsErrors=False,
                skipErrors=False,
                verbose=True
                )
        finally:
            log = open("map.log")
            self.info(log.read())
            log.close()

            self.info("--- Partial run time: %s" % self.get_partial_time())
            self.reset_partial_time()
            self.info("--- Comulative run time: %s\n\n" % self.get_total_time())

    def execute_makehess(self):
        """ execute makehess """
        self.info("\n\n\n= Running modelfit makehess =")

        cmd = (
            "ulimit -s unlimited && "
            #"python -u ${MWFITTING_DIR}/modelfit.py --config mwf.ini --makehess &> makehess.log"
            #"python -u /home/adriano.pieres/mwfitting/modelfit.py --config mwf.ini --makehess &> makehess.log"
            "python /home/adriano.pieres/mwf_git/mwfitting/multi_make_hess.py &> makehess.log"

        )

        try:
            osutil.OSUtils().run_eups_command(
                cmd,
                mw_fitting_multi_make_hess,
                stderrAsErrors=False,
                skipErrors=False,
                verbose=True
                )
        finally:
            log = open("makehess.log")
            self.info(log.read())
            log.close()

            self.info("--- Partial run time: %s" % self.get_partial_time())
            self.reset_partial_time()
            self.info("--- Comulative run time: %s\n\n" % self.get_total_time())

    def execute_emcee(self):
        """ execute emcee """
        self.info("\n\n\n= Running modelfit emcee =")

        cmd = (
            "ulimit -s unlimited && "
            #"python -u ${MWFITTING_DIR}/modelfit.py --config mwf.ini --emcee &> emcee.log"
            "python -u /home/adriano.pieres/mwf_git/mwfitting/modelfit.py --config mwf.ini --emcee &> emcee.log"
        )

        try:
            osutil.OSUtils().run_eups_command(
                cmd,
                mw_fitting_modelfit,
                stderrAsErrors=False,
                skipErrors=False,
                verbose=True
                )
        finally:
            log = open("emcee.log")
            self.info(log.read())
            log.close()

            self.info("--- Partial run time: %s" % self.get_partial_time())
            self.reset_partial_time()
            self.info("--- Comulative run time: %s\n\n" % self.get_total_time())

    def execute_multi_make_mock(self):
        """ execute multi_make_mock """
        self.info("\n\n\n= Running multi_make_mock =")

        cmd = (
            "ulimit -s unlimited && "
            #"python -u ${MWFITTING_DIR}/modelfit.py --config mwf.ini --emcee &> emcee.log"
            "python -u /home/adriano.pieres/mwf_git/mwfitting/multi_make_mock.py &> multi_make_mock.log"
        )

        try:
            osutil.OSUtils().run_eups_command(
                cmd,
                mw_fitting_modelfit,
                stderrAsErrors=False,
                skipErrors=False,
                verbose=True
                )
        finally:
            log = open("multi_make_mock.log")
            self.info(log.read())
            log.close()

            self.info("--- Partial run time: %s" % self.get_partial_time())
            self.reset_partial_time()
            self.info("--- Comulative run time: %s\n\n" % self.get_total_time())

    def execute_dist_stars(self):
        """ execute dist_stats """
        self.info("\n\n\n= Running dist_stars =")

        cmd = (
            "ulimit -s unlimited && "
            #"python3 -u ${MWFITTING_DIR}/dist_stars.py &> mock.log"
            "python3 -u /home/adriano.pieres/mwf_git/mwfitting/dist_stars.py &> dist_stars.log"
        )

        try:
            osutil.OSUtils().run_eups_command(
                cmd,
                mw_fitting_map,
                stderrAsErrors=False,
                skipErrors=False,
                verbose=True
                )
        finally:
            log = open("dist_stars.log")
            self.info(log.read())
            log.close()

            self.info("--- Partial run time: %s" % self.get_partial_time())
            self.reset_partial_time()
            self.info("--- Comulative run time: %s\n\n" % self.get_total_time())

    def create_product_log(self):
        """ create product log """
        self.info("\n\n\n= Creating product log =")
        conf = cpio.ComponentConfig()
        pio = cpio.ProvenanceIO()

        self.io.beginLog("Input parameters", id="0")
        self.io.beginSection()
        self.io.beginSubSection()
        
        columns = ['Configuration parameter', 'Description', 'Value', 'Unit']
        self.add_table_from_file('MWF_pars.dat', columns, ' ')
        
        data_file = open('MWF_pars.dat', 'r')
        data_lines = data_file.readlines()
        data_file.close()
        self.io.beginStats()

        self.info("Reading files for information")

        for i in range(len(data_lines)):
            line = data_lines[i]
            data = line.split()
            data[1].replace('_', ' ')
            if (i==0):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Nside for Hess Diagrams',_publish='True')
            elif (i==1):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Nside footprint map',_publish='True')
            elif (i==2):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Surveys in the comparison',_publish='True')
            elif (i==3):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Minimum fraction of the cell area',_publish='True')
            elif (i==4):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Brighter magnitude limit to each survey',_publish='True')
            elif (i==5):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Fainter magnitude limit to each survey',_publish='True')
            elif (i==6):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Magnitude step',_publish='True')
            elif (i==7):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Bluer color limit to each survey',_publish='True')
            elif (i==8):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Redder color limit to each survey',_publish='True')
            elif (i==9):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Color step',_publish='True')
            elif (i==10):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Comparison fields in each survey',_publish='True')
            elif (i==11):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Overfactor (how many times models are simulated)',_publish='True')
            elif (i==12):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Emcee walkers',_publish='True')
            elif (i==13):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Steps for walkers',_publish='True')
            elif (i==14):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Steps for the burning phase',_publish='True')
            elif (i==15):
                self.io.putParam(_value=data[2],_type='string',_id='code',_section='Input Parameters',_name='Temperature (step length)',_publish='True')

        self.io.endStats()
        self.io.endSubSection()
        self.io.endSection()
        self.io.endLog(id="0")

        self.io.beginLog("Pointings", id="1")
        self.io.beginSection()
        self.io.beginSubSection()

        self.io.beginTable(_class="paginator")

        data_file = open('fields_list.dat', 'r')
        data_lines = data_file.readlines()
        data_file.close()

        columns = [
            "(data - model)/data (%)",
            "Field (Survey+ID)",
            "RA (deg)",
            "DEC (deg)",
            "l (deg)",
            "b (deg)",
            "Star counts (data)",
            "Star counts (model)",
            "Total area (sq deg)",
            "Av (average)"
            ]

        for i, col in enumerate(columns):
            self.io.beginColumn(col)
            for j, line in enumerate(data_lines):
                self.io.beginRow(str(j))
                data = line.split()
                self.io.addStat("", data[i])
                self.io.endRow()
            self.io.endColumn()

        self.io.endTable()

        os.system('cp plots/*.png .')

        for png in glob.glob('*_initial_guess.png'):
            self.io.addPlot(png, 'Plots with initial parameters', '')
            png2 = png.replace('initial_guess', 'best_model')
            self.io.addPlot(png2, 'Plots with best model', '')

        self.io.endSubSection()
        '''
        self.io.beginSubSection()

        lhood_title = (
             'The on-sky distribution of the regions being fitted. ' +
             'An ID number is displayed along with the region\'s ' +
             'contribution to the total chi2. The size of each point indicates ' +
             'the corresponding area and the colorbar represents the total chi2.'
             )

        # os.system('cp plots/lhood_map.png lhood.png')

        # self.io.addPlot('lhood.png', 'Map of pointings on the sky', lhood_title)

        self.io.endSubSection()
        '''
        self.io.beginSubSection()

        for png in glob.glob('*stat.png'):
            pattern = re.compile('DES\d*')
            _id = pattern.search(png)
            _id = _id.group()

            self.io.addPlot(png, _id, '')

        self.io.endSubSection()

        self.io.endSection()
        self.io.endLog(id="1")

        self.io.beginLog("Coverage Map", id="2")
        self.io.beginSection()
        self.io.beginSubSection()
        self.io.addPlot('foot_plot_file_EQU.png', 'Title', 'caption')
        self.io.addPlot('foot_plot_file_GAL.png', 'Title', 'caption')
        self.io.endSubSection()
        self.io.endSection()
        self.io.endLog(id="2")

        self.io.beginLog("MCMC", id="3")
        self.io.beginSection()
        self.io.beginSubSection()
        os.system('cp plots/burner.png burner.png')
        self.io.addPlot('burner.png', 'Posterior Likelihoods for the fitting parameters', 'caption')
        self.io.endSubSection()
        self.io.beginSubSection()
        os.system('cp plots/line-time.png line-time.png')
        self.io.addPlot('line-time.png', 'Evolution of the walkers along the steps', 'caption')
        self.io.endSubSection()
        self.io.endSection()
        self.io.endLog(id="3")

        self.io.beginLog("Profiles", id="4")
        self.io.beginSection()
        self.io.beginSubSection()
        self.io.addPlot('halo_profile.png', '', 'caption')
        self.io.endSubSection()
        self.io.beginSubSection()
        self.io.addPlot('thickdisk_profile.png', '', 'caption')
        self.io.endSubSection()
        self.io.endSection()
        self.io.endLog(id="4")

        self.io.beginLog("Best points", id="5")
        self.io.beginSection()
        self.io.beginSubSection()

        columns = ['Parameter (name)', 'Description', 'Unit', 'Initial Guess', 'Best-Fitting', 'Error(+)', 'Error(-)',
                   'Best-Fitting/Guess', 'Error(+)/Guess', 'Error(-)/Guess', 'Prior(-)', 'Prior(+)', 'True value']
        self.add_table_from_file('best.dat', columns, ' ')

        d_file = open('best.dat', 'r')
        d_lines = d_file.readlines()
        d_file.close()

        for i in range(len(d_lines)):
            line = d_lines[i]
            data = line.split()
            data[1].replace('_', ' ')
            if (i==0):
                self.io.putParam(_value=data[3],_type='string',_id='code',_section='MW Input Parameters',_name='Halo density at Sun position (solar masses per cubic parsec)',_publish='True')
                self.io.putParam(_value=data[10],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Halo density (lower limit)',_publish='True')
                self.io.putParam(_value=data[11],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Halo density (upper limit)',_publish='True')
            if (i==1):
                self.io.putParam(_value=data[3],_type='string',_id='code',_section='MW Input Parameters',_name='Thick disk vertical scale (parsec)',_publish='True')
                self.io.putParam(_value=data[10],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Thick disk vertical scale (lower limit)',_publish='True')
                self.io.putParam(_value=data[11],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Thick disk vertical scale (upper limit)',_publish='True')
            if (i==2):
                self.io.putParam(_value=data[3],_type='string',_id='code',_section='MW Input Parameters',_name='Halo q',_publish='True')
                self.io.putParam(_value=data[10],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Halo q (lower limit)',_publish='True')
                self.io.putParam(_value=data[11],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Halo q (upper limit)',_publish='True')
            if (i==3):
                self.io.putParam(_value=data[3],_type='string',_id='code',_section='MW Input Parameters',_name='Halo exponent (n)',_publish='True')
                self.io.putParam(_value=data[10],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Halo n (lower limit)',_publish='True')
                self.io.putParam(_value=data[11],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Halo n (upper limit)',_publish='True')
            if (i==4):
                self.io.putParam(_value=data[3],_type='string',_id='code',_section='MW Input Parameters',_name='Thick disk scale radius (parsec)',_publish='True')
                self.io.putParam(_value=data[10],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Thick disk scale radius (lower limit)',_publish='True')
                self.io.putParam(_value=data[11],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Thick disk scale radius (upper limit)',_publish='True')
            if (i==5):
                self.io.putParam(_value=data[3],_type='string',_id='code',_section='MW Input Parameters',_name='Thick disk surface density at Sun position (solar masses per square parsec)',_publish='True')
                self.io.putParam(_value=data[10],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Thick disk density (lower limit)',_publish='True')
                self.io.putParam(_value=data[11],_type='string',_id='code',_section='MW Parameters Priors (relative to input parameters)',_name='Thick disk density (upper limit)',_publish='True')

        self.io.endSubSection()

        self.io.endSection()
        self.io.endLog(id="5")


        self.io.beginLog("Density maps for simulated catalog", id="6")
        self.io.beginSection()

        if os.path.isfile('plots/des_diff_map_total_EQU.png'):
            os.system('cp plots/des_diff_map_total_EQU.png des_diff_map_total_EQU.png')
            os.system('cp plots/des_diff_map_total_GAL.png des_diff_map_total_GAL.png')
            os.system('cp plots/des_diff_map_0_EQU.png des_diff_map_0_EQU.png')
            os.system('cp plots/des_diff_map_1_EQU.png des_diff_map_1_EQU.png')
            os.system('cp plots/des_diff_map_2_EQU.png des_diff_map_2_EQU.png')
            os.system('cp plots/des_diff_map_3_EQU.png des_diff_map_3_EQU.png')
            os.system('cp plots/des_diff_map_0_GAL.png des_diff_map_0_GAL.png')
            os.system('cp plots/des_diff_map_1_GAL.png des_diff_map_1_GAL.png')
            os.system('cp plots/des_diff_map_2_GAL.png des_diff_map_2_GAL.png')
            os.system('cp plots/des_diff_map_3_GAL.png des_diff_map_3_GAL.png')

            self.info("\n\n\n= Showing density maps for DES =")
 
            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_total_EQU.png', 'Density maps for full color range', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_total_GAL.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_0_EQU.png', 'Density maps for specific color ranges', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_1_EQU.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_2_EQU.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_3_EQU.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_0_GAL.png', 'Density maps for specific color ranges', 'caption')
            self.io.endSubSection()

            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_1_GAL.png', '', 'caption')
            self.io.endSubSection()

            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_2_GAL.png', '', 'caption')
            self.io.endSubSection()

            self.io.beginSubSection()
            self.io.addPlot('des_diff_map_3_GAL.png', '', 'caption')
            self.io.endSubSection()
        else:
            self.info("\n\n\n= No density maps for DES =")
            pass

        if os.path.isfile('plots/sdss_diff_map_total_EQU.png'):
            os.system('cp plots/sdss_diff_map_total_EQU.png sdss_diff_map_total_EQU.png')
            os.system('cp plots/sdss_diff_map_total_GAL.png sdss_diff_map_total_GAL.png')
            os.system('cp plots/sdss_diff_map_0_EQU.png sdss_diff_map_0_EQU.png')
            os.system('cp plots/sdss_diff_map_1_EQU.png sdss_diff_map_1_EQU.png')
            os.system('cp plots/sdss_diff_map_2_EQU.png sdss_diff_map_2_EQU.png')
            os.system('cp plots/sdss_diff_map_3_EQU.png sdss_diff_map_3_EQU.png')
            os.system('cp plots/sdss_diff_map_0_GAL.png sdss_diff_map_0_GAL.png')
            os.system('cp plots/sdss_diff_map_1_GAL.png sdss_diff_map_1_GAL.png')
            os.system('cp plots/sdss_diff_map_2_GAL.png sdss_diff_map_2_GAL.png')
            os.system('cp plots/sdss_diff_map_3_GAL.png sdss_diff_map_3_GAL.png')

            self.info("\n\n\n= Showing density maps for SDSS =")
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_total_EQU.png', 'Density maps for full color range', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_total_GAL.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_0_EQU.png', 'Density maps for specific color ranges', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_1_EQU.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_2_EQU.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_3_EQU.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_0_GAL.png', 'Density maps for specific color ranges', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_1_GAL.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_2_GAL.png', '', 'caption')
            self.io.endSubSection()
 
            self.io.beginSubSection()
            self.io.addPlot('sdss_diff_map_3_GAL.png', '', 'caption')
            self.io.endSubSection()
        else:
            self.info("\n\n\n= No density maps for SDSS =")
            pass
     
        products = pio.products
        curr_direc = os.getcwd()

        if os.path.isfile(curr_direc + '/data/des_mockcat_pos_mag.fits'):
            self.info('\n\n\nAFTER= %s' % curr_direc)
            MWF_SIM_VAC = self.io.putFile(_value=curr_direc + '/data/des_mockcat_pos_mag.fits',
                         _class='vac_simulated_mwfitting',
                         _mimetype='application/fits-table',
                         _publish='True')
        else:
            self.info("\n\n\n= No simulated catalog maps for DES =")
            pass
        if os.path.isfile(curr_direc + '/data/sdss_mockcat_pos_mag.fits'):
            self.info('\n\n\n= %s' % curr_direc)
            MWF_SIM_VAC = self.io.putFile(_value=curr_direc + '/data/sdss_mockcat_pos_mag.fits',
                         _class='vac_simulated_mwfitting',
                         _mimetype='application/fits-table',
                         _publish='True')
        else:
            self.info("\n\n\n= No simulated catalog maps for SDSS =")
            pass

        self.io.endSection()
        self.io.endLog(id="6")

def run(debug=True):
    component = MwFitting(orchestration.io.ComponentIO(), debug)
    component.start()


