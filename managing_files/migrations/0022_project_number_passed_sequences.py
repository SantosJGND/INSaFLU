# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2021-04-21 17:58
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('managing_files', '0021_projectsample_seq_name_all_consensus'),
    ]

    operations = [
        migrations.AddField(
            model_name='project',
            name='number_passed_sequences',
            field=models.SmallIntegerField(default=-1),
        ),
    ]
