# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2020-10-18 10:36
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('managing_files', '0015_auto_20200323_2050'),
    ]

    operations = [
        migrations.AddField(
            model_name='projectsample',
            name='is_mask_consensus_sequences',
            field=models.BooleanField(default=False),
        ),
    ]
