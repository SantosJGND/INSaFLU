"""
Django settings for fluwebvirus project.

Generated by 'django-admin startproject' using Django 1.11.6.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.11/ref/settings/
"""

import os
from decouple import config

### running tests in command line
RUN_TEST_IN_COMMAND_LINE = True	## used in abricate

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = config('SECRET_KEY')

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = config('DEBUG', default=True, cast=bool)

ALLOWED_HOSTS = config('ALLOWED_HOSTS', cast=lambda v: [s.strip() for s in v.split(',')])

## google recaptcha
GOOGLE_RECAPTCHA_SECRET_KEY = config('GOOGLE_RECAPTCHA_SECRET_KEY')
SITE_KEY = config('SITE_KEY')

### crispy template
CRISPY_TEMPLATE_PACK = 'bootstrap4'
BREADCRUMBS_TEMPLATE = "django_bootstrap_breadcrumbs/bootstrap4.html"

CSRF_COOKIE_AGE = None
CSRF_COOKIE_DOMAIN = '.min-saude.pt'
CSRF_COOKIE_HTTPONLY = True
CSRF_COOKIE_SECURE = True
CSRF_USE_SESSIONS = True

### threads to use in several software
THREADS_TO_RUN_FAST= config('THREADS_TO_RUN_FAST', default=3, cast=int)
THREADS_TO_RUN_SLOW = config('THREADS_TO_RUN_SLOW', default=3, cast=int)

#https://www.digitalocean.com/community/tutorials/how-to-create-an-ssl-certificate-on-apache-for-centos-7
#https://gist.github.com/bradmontgomery/6487319
#SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')
if (config('SECURE_SSL_REDIRECT', default=False, cast=bool)):
	SECURE_SSL_REDIRECT = True
	SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')

## add google analytics
ADD_GOOGLE_ANALYTICS = config('ADD_GOOGLE_ANALYTICS', default=False, cast=bool)

## to show login anonymous
SHOW_LOGIN_ANONYMOUS = config('SHOW_LOGIN_ANONYMOUS', default=False, cast=bool)

## make the doen size of the fastq files to 50MB
## if the DOWN_SIZE_FASTQ_FILES is false the maximum fastq input files is 50MB
DOWN_SIZE_FASTQ_FILES = config('DOWN_SIZE_FASTQ_FILES', default=False, cast=bool)

## run process in SGE, otherwise run in qcluster
RUN_SGE  = config('RUN_SGE', default=False, cast=bool)
SGE_ROOT = config('SGE_ROOT')

# Application definition
INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.gis',
    'crispy_forms',
    'crispy_forms_foundation',
    'django_tables2',
    'bootstrap4',
    'django_q',
    'django_user_agents',
    'django_bootstrap_breadcrumbs',
    'managing_files.apps.ManagingFilesConfig',
    'manage_virus.apps.ManageVirusConfig',
    'phylogeny.apps.PhylogenyConfig',
    'log_login.apps.LogLoginConfig',
    'extend_user.apps.ExtendUserConfig',
    'crequest',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django_user_agents.middleware.UserAgentMiddleware',
    'crequest.middleware.CrequestMiddleware',
]

ROOT_URLCONF = 'fluwebvirus.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, 'templates'),
            ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

## $ python3 manage.py qcluster
## $ python3 manage.py qmonitor
## $ python3 manage.py qinfo
## this configuraton doesn't work
# Q_CLUSTER = {
#     'name': 'insaFlu',
#     'workers': 1,	## number of queues, some problems with database, need to bee always at one
#     'recycle': 500,
#     'compress': True,
#     'cached': False,
#     'save_limit': 250,
#     'queue_limit': 1,
#     'cpu_affinity': 2,	## number of processors by queue
#     'catch_up': False,	# Ignore un-run scheduled tasks
#     'label': 'Django Q',
#     'orm': 'default'
# }

#redis defaults
Q_CLUSTER = {
	'name': 'insaFlu',
	'workers': 1,	## number of queues, some problems with database, need to bee always at one
	'recycle': 500,
	'compress': True,
	'save_limit': 250,
    'cpu_affinity': 2,	## number of processors by queue
    'catch_up': False,	# Ignore un-run scheduled tasks
    'label': 'Django Q',

    'redis': {
        'host': 'localhost',
        'port': 6379,
        'db': 0,
        'password': None,
        'socket_timeout': None,
        'charset': 'utf-8',
        'errors': 'strict',
        'unix_socket_path': None
    }
}

CACHES = {
    'default': {
		'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'djangoq-localmem',
    }
}

#### EMAIL
#EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
## NTLM in  preferred_auths = [AUTH_CRAM_MD5, AUTH_PLAIN, AUTH_LOGIN] 
# sudo ssh -p 2023 -L 25:192.168.32.99:25 insa@193.137.95.75

# http://www.techspacekh.com/configuring-postfix-to-relay-mail-to-local-exchange-mail-server-in-rhel-centos-7/
# http://www.postfix.org/SASL_README.html#saslauthd
## see more dropbox/insa
EMAIL_BACKEND = config('EMAIL_BACKEND')
EMAIL_HOST = '127.0.0.1' ###config('EMAIL_HOST')
DEFAULT_FROM_EMAIL = config('EMAIL_NAME')
EMAIL_PORT = config('EMAIL_PORT', cast=int)
EMAIL_HOST_USER = ''   		### config('EMAIL_HOST_USER')
EMAIL_HOST_PASSWORD = '' 	### config('EMAIL_HOST_PASSWORD')
EMAIL_USE_TLS = config('EMAIL_USE_TLS', cast=bool)

WSGI_APPLICATION = 'fluwebvirus.wsgi.application'

# Name of cache backend to cache user agents. If it not specified default
# cache alias will be used. Set to `None` to disable caching.
USER_AGENTS_CACHE = 'default'

# Database
# https://docs.djangoproject.com/en/1.11/ref/settings/#databases
# https://www.digitalocean.com/community/tutorials/how-to-use-postgresql-with-your-django-application-on-ubuntu-14-04
# DATABASES = {
#     'default': {
#         'ENGINE': 'django.db.backends.sqlite3',
#         'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
#     }
# }


## to reuse DB 
# os.environ['REUSE_DB'] = "1"
DATABASES = {
    'default': {
        'ENGINE': 'django.contrib.gis.db.backends.postgis',
##        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'fluwebvirus_test',
        'USER': config('DB_USER'),
        'PASSWORD': config('DB_PASSWORD'),
        'HOST': config('DB_HOST'),
        'PORT': '5432',
    },
}

# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]



# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/
#DATE_INPUT_FORMATS = ('%d-%m-%Y','%Y-%m-%d')

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Lisbon'

DATE_FORMAT_FOR_TABLE = '%d-%m-%Y %H:%M'
DATETIME_FORMAT = '%d-%m-%Y %H:%M'
DATETIME_INPUT_FORMATS = ['%d-%m-%Y', '%d/%m/%Y']	## it's necessary to look which kind of date is returned from forms to correct the format

## set the format in forms.DateField
## This will only work if USE_L10N is False. You may also need to set DATE_FORMAT used when printing a date in the templates
USE_L10N = False

USE_TZ = False	## enable time zone



# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.11/howto/static-files/

### Limit the siz files
# Look at the LimitRequestBody directive.
# https://stackoverflow.com/questions/2472422/django-file-upload-size-limit
# https://stackoverflow.com/questions/5871730/need-a-minimal-django-file-upload-example

# STATICFILES_DIRS is the list of folder where Django will search for additional static files, in addition to each static folder of each app installed.
STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "static"),
)
# STATIC_ROOT is the folder where every static files will be stored after a manage.py collectstatic
STATIC_ROOT = os.path.join(BASE_DIR, 'static_all')	## is the absolute path to the directory where collectstatic will collect static files for deployment.
STATIC_URL = '/static/' 		## is the URL to use when referring to static files located in STATIC_ROOT.

MEDIA_ROOT_TEST = os.path.join('/tmp/tests_insa_flu')
MEDIA_ROOT = os.path.join(BASE_DIR, '/tmp/tests_insa_flu')
MEDIA_URL = '/media/'


# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'django.server': {
            '()': 'django.utils.log.ServerFormatter',
            'format': '[%(server_time)s] %(message)s',
        },
		'verbose': {
        	'format': '%(levelname)s %(asctime)s %(module)s: %(message)s'
    	}
    },
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        },
        'require_debug_true': {
            '()': 'django.utils.log.RequireDebugTrue',
        },
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'class': 'django.utils.log.AdminEmailHandler'
        },
		'file_debug': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': '/var/log/insaFlu/debug.log',
            'formatter': 'verbose',
        },
		'file_warning': {
            'level': 'WARNING',
            'class': 'logging.FileHandler',
            'filename': '/var/log/insaFlu/warning.log',
            'formatter': 'verbose',
        },
		'console':{
            'level':'DEBUG',
            'class':'logging.StreamHandler',
            'formatter': 'verbose'
        },
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'filters': ['require_debug_false'],
            'level': 'ERROR',
            'propagate': True,
        },
		'fluWebVirus.debug': {
            'handlers': ['file_debug'],
            'filters': ['require_debug_true'],
            'level': 'DEBUG',
            'propagate': True,
        },
		'fluWebVirus.production': {
            'handlers': ['file_warning'],
            'filters': ['require_debug_false'],
            'level': 'WARNING',		## third level of log
            'propagate': True,
        },
    }
}